package uk.ac.manchester.tornado.examples.compute;

import java.util.Arrays;
import java.util.stream.IntStream;

import uk.ac.manchester.tornado.api.GridScheduler;
import uk.ac.manchester.tornado.api.ImmutableTaskGraph;
import uk.ac.manchester.tornado.api.KernelContext;
import uk.ac.manchester.tornado.api.TaskGraph;
import uk.ac.manchester.tornado.api.TornadoExecutionPlan;
import uk.ac.manchester.tornado.api.WorkerGrid1D;
import uk.ac.manchester.tornado.api.enums.DataTransferMode;
import uk.ac.manchester.tornado.api.math.TornadoMath;
import uk.ac.manchester.tornado.api.types.arrays.FloatArray;

public class RMSSIMDReduction {
    private static final int   LOCAL_SIZE = 256;   // 8 warps per workgroup
    private static final float ERMS_NORM  = 1e-5f;
    private static final int   WARMUP     = 50;
    private static final int   ITERATIONS = 200;

    // =========================================================================
    // Kernel 1: original — threadgroup memory tree reduction
    // =========================================================================
    private static void reductionOriginal(
            KernelContext context, FloatArray output, FloatArray x,
            int size, float ermsNorm, int localMemSize) {

        int gid       = context.globalIdx;
        int lid       = context.localIdx;
        int groupId   = context.groupIdx;
        int groupSize = context.localGroupSizeX;

        float[] localX = context.allocateFloatLocalArray(localMemSize);

        if (gid < size) {
            localX[lid] = x.get(gid);
            localX[lid] = localX[lid] * localX[lid];
        } else {
            localX[lid] = 0.0f;
        }

        for (int stride = groupSize / 2; stride > 0; stride /= 2) {
            context.localBarrier();
            if (lid < stride) {
                localX[lid] += localX[lid + stride];
            }
        }

        if (lid == 0) {
            output.set(groupId + 1, localX[0]);
        }

        if (gid == 0) {
            float ss = 0.0f;
            for (int i = 1; i <= (size / localMemSize); i++) {
                ss += output.get(i);
            }
            ss /= size;
            ss += ermsNorm;
            ss = 1.0f / TornadoMath.sqrt(ss);
            output.set(0, ss);
        }
    }

    // =========================================================================
    // Kernel 2: Shuffle + Shuffle
    // Level 1: manual simdShuffleDown butterfly (= warp_reduce_sum in CUDA ref)
    // Level 2: manual simdShuffleDown butterfly on the per-warp sums
    // =========================================================================
    private static void reductionShuffleShuffle(
            KernelContext context, FloatArray output, FloatArray x,
            int size, float ermsNorm, int localMemSize) {

        int gid       = context.globalIdx;
        int lid       = context.localIdx;
        int groupId   = context.groupIdx;
        int groupSize = context.localGroupSizeX;

        // Level 1: warp butterfly
        float val = (gid < size) ? x.get(gid) : 0.0f;
        val = val * val;
        val += context.simdShuffleDown(val, 16);
        val += context.simdShuffleDown(val, 8);
        val += context.simdShuffleDown(val, 4);
        val += context.simdShuffleDown(val, 2);
        val += context.simdShuffleDown(val, 1);

        float[] sharedSums = context.allocateFloatLocalArray(32);
        int lane   = lid % 32;
        int warpId = lid / 32;
        if (lane == 0) sharedSums[warpId] = val;
        context.localBarrier();

        // Level 2: first warp reduces the warp sums via butterfly
        int numWarps = groupSize / 32;
        float warpSum = (lid < numWarps) ? sharedSums[lane] : 0.0f;
        warpSum += context.simdShuffleDown(warpSum, 16);
        warpSum += context.simdShuffleDown(warpSum, 8);
        warpSum += context.simdShuffleDown(warpSum, 4);
        warpSum += context.simdShuffleDown(warpSum, 2);
        warpSum += context.simdShuffleDown(warpSum, 1);

        if (lid == 0) output.set(groupId + 1, warpSum);

        if (gid == 0) {
            float ss = 0.0f;
            for (int i = 1; i <= (size / localMemSize); i++) {
                ss += output.get(i);
            }
            ss /= size;
            ss += ermsNorm;
            ss = 1.0f / TornadoMath.sqrt(ss);
            output.set(0, ss);
        }
    }

    // =========================================================================
    // Kernel 3: simdSum + Shuffle
    // Level 1: simdSum (compiler-owned — identical output to K2 level 1)
    // Level 2: manual simdShuffleDown butterfly on the per-warp sums
    // =========================================================================
    private static void reductionSimdSumShuffle(
            KernelContext context, FloatArray output, FloatArray x,
            int size, float ermsNorm, int localMemSize) {

        int gid       = context.globalIdx;
        int lid       = context.localIdx;
        int groupId   = context.groupIdx;
        int groupSize = context.localGroupSizeX;

        // Level 1: simdSum per warp
        float val = (gid < size) ? x.get(gid) : 0.0f;
        val = val * val;
        float subgroupSum = context.simdSum(val);

        float[] sharedSums = context.allocateFloatLocalArray(32);
        int lane   = lid % 32;
        int warpId = lid / 32;
        if (lane == 0) sharedSums[warpId] = subgroupSum;
        context.localBarrier();

        // Level 2: manual butterfly on the warp sums
        int numWarps = groupSize / 32;
        float warpSum = (lid < numWarps) ? sharedSums[lane] : 0.0f;
        warpSum += context.simdShuffleDown(warpSum, 16);
        warpSum += context.simdShuffleDown(warpSum, 8);
        warpSum += context.simdShuffleDown(warpSum, 4);
        warpSum += context.simdShuffleDown(warpSum, 2);
        warpSum += context.simdShuffleDown(warpSum, 1);

        if (lid == 0) output.set(groupId + 1, warpSum);

        if (gid == 0) {
            float ss = 0.0f;
            for (int i = 1; i <= (size / localMemSize); i++) {
                ss += output.get(i);
            }
            ss /= size;
            ss += ermsNorm;
            ss = 1.0f / TornadoMath.sqrt(ss);
            output.set(0, ss);
        }
    }

    // =========================================================================
    // Kernel 4: simdSum + simdSum
    // Level 1: simdSum per warp
    // Level 2: simdSum on the numWarps values loaded into the first warp
    // =========================================================================
    private static void reductionSimdSumSimdSum(
            KernelContext context, FloatArray output, FloatArray x,
            int size, float ermsNorm, int localMemSize) {

        int gid       = context.globalIdx;
        int lid       = context.localIdx;
        int groupId   = context.groupIdx;
        int groupSize = context.localGroupSizeX;

        // Level 1: simdSum per warp
        float val = (gid < size) ? x.get(gid) : 0.0f;
        val = val * val;
        float subgroupSum = context.simdSum(val);

        float[] sharedSums = context.allocateFloatLocalArray(32);
        int lane   = lid % 32;
        int warpId = lid / 32;
        if (lane == 0) sharedSums[warpId] = subgroupSum;
        context.localBarrier();

        // Level 2: simdSum on the warp sums (lanes beyond numWarps get 0.0f)
        int numWarps = groupSize / 32;
        float warpVal    = (lid < numWarps) ? sharedSums[lane] : 0.0f;
        float totalSum   = context.simdSum(warpVal);

        if (lid == 0) output.set(groupId + 1, totalSum);

        if (gid == 0) {
            float ss = 0.0f;
            for (int i = 1; i <= (size / localMemSize); i++) {
                ss += output.get(i);
            }
            ss /= size;
            ss += ermsNorm;
            ss = 1.0f / TornadoMath.sqrt(ss);
            output.set(0, ss);
        }
    }

    // =========================================================================
    // Sequential reference
    // =========================================================================
    private static float sequentialRmsNormFactor(FloatArray x, int size, float ermsNorm) {
        double ss = 0.0;
        for (int i = 0; i < size; i++) {
            float v = x.get(i);
            ss += (double) v * v;
        }
        ss /= size;
        ss += ermsNorm;
        return (float) (1.0 / Math.sqrt(ss));
    }

    // =========================================================================
    // Helpers
    // =========================================================================
    private static long[] benchmark(TornadoExecutionPlan plan, GridScheduler grid) {
        for (int i = 0; i < WARMUP; i++) plan.withGridScheduler(grid).execute();
        long[] times = new long[ITERATIONS];
        for (int i = 0; i < ITERATIONS; i++) {
            long t0 = System.nanoTime();
            plan.withGridScheduler(grid).execute();
            times[i] = System.nanoTime() - t0;
        }
        return times;
    }

    private static void printStats(String label, long[] timesNs) {
        double[] ms = Arrays.stream(timesNs).mapToDouble(t -> t / 1e6).toArray();
        double avg = Arrays.stream(ms).average().orElse(0);
        double min = Arrays.stream(ms).min().orElse(0);
        double max = Arrays.stream(ms).max().orElse(0);
        double std = Math.sqrt(Arrays.stream(ms)
                .map(v -> (v - avg) * (v - avg)).average().orElse(0));
        System.out.printf("  %-38s  avg=%6.3f ms  min=%6.3f ms  max=%6.3f ms  std=%5.3f ms%n",
                label, avg, min, max, std);
    }

    private static void checkResult(String label, float expected, float actual) {
        float relErr = Math.abs(actual - expected) / Math.abs(expected);
        System.out.printf("  %-38s  result=%.8f  relErr=%.2e  %s%n",
                label, actual, relErr, relErr < 1e-3f ? "✓" : "✗ MISMATCH");
    }

    private static String resultMark(float expected, float actual) {
        float relErr = (expected == 0) ? Math.abs(actual)
                : Math.abs(actual - expected) / Math.abs(expected);
        return relErr < 1e-3f
                ? String.format("%.6f ✓", actual)
                : String.format("%.6f ✗ (exp %.6f)", actual, expected);
    }

    // =========================================================================
    // main
    // =========================================================================
    public static void main(String[] args) throws Exception {

        final int size      = 2048;
        final int numGroups = size / LOCAL_SIZE;

        System.out.println("RMS-norm Reduction — Four-Way Strategy Comparison");
        System.out.println("===================================================");
        System.out.printf("dim=%d  localSize=%d (%d warps)  groups=%d  ermsNorm=%e%n%n",
                size, LOCAL_SIZE, LOCAL_SIZE / 32, numGroups, ERMS_NORM);

        FloatArray input = new FloatArray(size);
        IntStream.range(0, size).forEach(i -> input.set(i, (float) Math.sin(i)));
        float expected = sequentialRmsNormFactor(input, size, ERMS_NORM);
        System.out.printf("Sequential reference = %.8f%n%n", expected);

        FloatArray outOrig = new FloatArray(numGroups + 1);
        FloatArray outShSh = new FloatArray(numGroups + 1);
        FloatArray outSsSh = new FloatArray(numGroups + 1);
        FloatArray outSsSs = new FloatArray(numGroups + 1);

        KernelContext ctx = new KernelContext();

        // ---- Kernel 1: threadgroup tree ----
        WorkerGrid1D gridOrig = new WorkerGrid1D(size);
        gridOrig.setLocalWork(LOCAL_SIZE, 1, 1);
        GridScheduler schedulerOrig = new GridScheduler("orig.t0", gridOrig);

        TaskGraph origGraph = new TaskGraph("orig")
                .transferToDevice(DataTransferMode.FIRST_EXECUTION, input)
                .task("t0", RMSSIMDReduction::reductionOriginal,
                        ctx, outOrig, input, size, ERMS_NORM, LOCAL_SIZE)
                .transferToHost(DataTransferMode.EVERY_EXECUTION, outOrig);
        ImmutableTaskGraph origITG = origGraph.snapshot();

        // ---- Kernel 2: Shuffle + Shuffle ----
        WorkerGrid1D gridShSh = new WorkerGrid1D(size);
        gridShSh.setLocalWork(LOCAL_SIZE, 1, 1);
        GridScheduler schedulerShSh = new GridScheduler("shsh.t0", gridShSh);

        TaskGraph shShGraph = new TaskGraph("shsh")
                .transferToDevice(DataTransferMode.FIRST_EXECUTION, input)
                .task("t0", RMSSIMDReduction::reductionShuffleShuffle,
                        ctx, outShSh, input, size, ERMS_NORM, LOCAL_SIZE)
                .transferToHost(DataTransferMode.EVERY_EXECUTION, outShSh);
        ImmutableTaskGraph shShITG = shShGraph.snapshot();

        // ---- Kernel 3: simdSum + Shuffle ----
        WorkerGrid1D gridSsSh = new WorkerGrid1D(size);
        gridSsSh.setLocalWork(LOCAL_SIZE, 1, 1);
        GridScheduler schedulerSsSh = new GridScheduler("sssh.t0", gridSsSh);

        TaskGraph ssShGraph = new TaskGraph("sssh")
                .transferToDevice(DataTransferMode.FIRST_EXECUTION, input)
                .task("t0", RMSSIMDReduction::reductionSimdSumShuffle,
                        ctx, outSsSh, input, size, ERMS_NORM, LOCAL_SIZE)
                .transferToHost(DataTransferMode.EVERY_EXECUTION, outSsSh);
        ImmutableTaskGraph ssShITG = ssShGraph.snapshot();

        // ---- Kernel 4: simdSum + simdSum ----
        WorkerGrid1D gridSsSs = new WorkerGrid1D(size);
        gridSsSs.setLocalWork(LOCAL_SIZE, 1, 1);
        GridScheduler schedulerSsSs = new GridScheduler("ssss.t0", gridSsSs);

        TaskGraph ssSsGraph = new TaskGraph("ssss")
                .transferToDevice(DataTransferMode.FIRST_EXECUTION, input)
                .task("t0", RMSSIMDReduction::reductionSimdSumSimdSum,
                        ctx, outSsSs, input, size, ERMS_NORM, LOCAL_SIZE)
                .transferToHost(DataTransferMode.EVERY_EXECUTION, outSsSs);
        ImmutableTaskGraph ssSsITG = ssSsGraph.snapshot();

        // ---- Run benchmarks ----
        long[] tOrig, tShSh, tSsSh, tSsSs;
        try (TornadoExecutionPlan planOrig = new TornadoExecutionPlan(origITG);
             TornadoExecutionPlan planShSh = new TornadoExecutionPlan(shShITG);
             TornadoExecutionPlan planSsSh = new TornadoExecutionPlan(ssShITG);
             TornadoExecutionPlan planSsSs = new TornadoExecutionPlan(ssSsITG)) {

            System.out.println("Running benchmarks...");
            tOrig = benchmark(planOrig, schedulerOrig);
            tShSh = benchmark(planShSh, schedulerShSh);
            tSsSh = benchmark(planSsSh, schedulerSsSh);
            tSsSs = benchmark(planSsSs, schedulerSsSs);
        }

        // ---- Correctness ----
        System.out.println("\nCorrectness");
        System.out.println("-----------");
        checkResult("1. Threadgroup tree (original)",  expected, outOrig.get(0));
        checkResult("2. Shuffle + Shuffle (CUDA ref)", expected, outShSh.get(0));
        checkResult("3. simdSum + Shuffle",            expected, outSsSh.get(0));
        checkResult("4. simdSum + simdSum",            expected, outSsSs.get(0));

        // ---- Performance ----
        System.out.println("\nPerformance (end-to-end dispatch + kernel + readback)");
        System.out.println("------------------------------------------------------");
        printStats("1. Threadgroup tree (original)",  tOrig);
        printStats("2. Shuffle + Shuffle (CUDA ref)", tShSh);
        printStats("3. simdSum + Shuffle",            tSsSh);
        printStats("4. simdSum + simdSum",            tSsSs);

        double avgOrig = Arrays.stream(tOrig).average().orElse(1);
        double avgShSh = Arrays.stream(tShSh).average().orElse(1);
        double avgSsSh = Arrays.stream(tSsSh).average().orElse(1);
        double avgSsSs = Arrays.stream(tSsSs).average().orElse(1);

        System.out.println("\nSpeedup vs threadgroup tree:");
        System.out.printf("  2. Shuffle + Shuffle : %.2fx%n", avgOrig / avgShSh);
        System.out.printf("  3. simdSum + Shuffle : %.2fx%n", avgOrig / avgSsSh);
        System.out.printf("  4. simdSum + simdSum : %.2fx%n", avgOrig / avgSsSs);

        double best = Math.min(avgShSh, Math.min(avgSsSh, avgSsSs));
        String bestName = best == avgShSh ? "2. Shuffle + Shuffle"
                : best == avgSsSh ? "3. simdSum + Shuffle"
                :                   "4. simdSum + simdSum";
        System.out.printf("%nFastest: %s%n", bestName);

        // ---- Correctness sweep over realistic model dims ----
        System.out.println("\nModel dimension correctness sweep");
        System.out.println("---------------------------------");
        System.out.printf("  %-20s  %-8s  %-20s  %-20s  %-20s  %-20s%n",
                "model/dim", "groups", "1.tree", "2.shuf+shuf", "3.sum+shuf", "4.sum+sum");

        int[]    dims   = {512, 1024, 2048, 3072, 3584, 4096};
        String[] models = {"test/512", "test/1024",
                "Llama3.2-1B/2048", "Llama3.2-3B/3072",
                "Qwen2.5-7B/3584",  "Llama3.1-8B/4096"};

        KernelContext ctx2 = new KernelContext();

        for (int m = 0; m < dims.length; m++) {
            int dim    = dims[m];
            int padded = ((dim + LOCAL_SIZE - 1) / LOCAL_SIZE) * LOCAL_SIZE;
            int ng     = padded / LOCAL_SIZE;

            FloatArray in2 = new FloatArray(padded);
            IntStream.range(0, dim).forEach(i -> in2.set(i, (float) Math.sin(i)));
            float exp2 = sequentialRmsNormFactor(in2, padded, ERMS_NORM);

            FloatArray o1 = new FloatArray(ng + 1);
            FloatArray o2 = new FloatArray(ng + 1);
            FloatArray o3 = new FloatArray(ng + 1);
            FloatArray o4 = new FloatArray(ng + 1);

            String pfx = "d" + m;

            WorkerGrid1D w1 = new WorkerGrid1D(padded); w1.setLocalWork(LOCAL_SIZE, 1, 1);
            WorkerGrid1D w2 = new WorkerGrid1D(padded); w2.setLocalWork(LOCAL_SIZE, 1, 1);
            WorkerGrid1D w3 = new WorkerGrid1D(padded); w3.setLocalWork(LOCAL_SIZE, 1, 1);
            WorkerGrid1D w4 = new WorkerGrid1D(padded); w4.setLocalWork(LOCAL_SIZE, 1, 1);

            GridScheduler s1 = new GridScheduler(pfx + "orig.t0", w1);
            GridScheduler s2 = new GridScheduler(pfx + "shsh.t0", w2);
            GridScheduler s3 = new GridScheduler(pfx + "sssh.t0", w3);
            GridScheduler s4 = new GridScheduler(pfx + "ssss.t0", w4);

            try (TornadoExecutionPlan p1 = new TornadoExecutionPlan(new TaskGraph(pfx + "orig")
                    .transferToDevice(DataTransferMode.EVERY_EXECUTION, in2)
                    .task("t0", RMSSIMDReduction::reductionOriginal,
                            ctx2, o1, in2, padded, ERMS_NORM, LOCAL_SIZE)
                    .transferToHost(DataTransferMode.EVERY_EXECUTION, o1).snapshot());
                 TornadoExecutionPlan p2 = new TornadoExecutionPlan(new TaskGraph(pfx + "shsh")
                         .transferToDevice(DataTransferMode.EVERY_EXECUTION, in2)
                         .task("t0", RMSSIMDReduction::reductionShuffleShuffle,
                                 ctx2, o2, in2, padded, ERMS_NORM, LOCAL_SIZE)
                         .transferToHost(DataTransferMode.EVERY_EXECUTION, o2).snapshot());
                 TornadoExecutionPlan p3 = new TornadoExecutionPlan(new TaskGraph(pfx + "sssh")
                         .transferToDevice(DataTransferMode.EVERY_EXECUTION, in2)
                         .task("t0", RMSSIMDReduction::reductionSimdSumShuffle,
                                 ctx2, o3, in2, padded, ERMS_NORM, LOCAL_SIZE)
                         .transferToHost(DataTransferMode.EVERY_EXECUTION, o3).snapshot());
                 TornadoExecutionPlan p4 = new TornadoExecutionPlan(new TaskGraph(pfx + "ssss")
                         .transferToDevice(DataTransferMode.EVERY_EXECUTION, in2)
                         .task("t0", RMSSIMDReduction::reductionSimdSumSimdSum,
                                 ctx2, o4, in2, padded, ERMS_NORM, LOCAL_SIZE)
                         .transferToHost(DataTransferMode.EVERY_EXECUTION, o4).snapshot())) {
                p1.withGridScheduler(s1).execute();
                p2.withGridScheduler(s2).execute();
                p3.withGridScheduler(s3).execute();
                p4.withGridScheduler(s4).execute();
            }

            System.out.printf("  %-20s  %-8d  %-20s  %-20s  %-20s  %-20s%n",
                    models[m], ng,
                    resultMark(exp2, o1.get(0)), resultMark(exp2, o2.get(0)),
                    resultMark(exp2, o3.get(0)), resultMark(exp2, o4.get(0)));
        }
    }
}
