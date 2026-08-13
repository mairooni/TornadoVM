package perfprobe;

import java.util.ArrayList;
import java.util.List;

import uk.ac.manchester.tornado.api.ImmutableTaskGraph;
import uk.ac.manchester.tornado.api.TaskGraph;
import uk.ac.manchester.tornado.api.TornadoExecutionPlan;
import uk.ac.manchester.tornado.api.enums.DataTransferMode;
import uk.ac.manchester.tornado.api.types.arrays.FloatArray;

/**
 * Prices TornadoVM's kernel-compilation path across four kernel complexity tiers, so the cost of
 * the JVMCI removal can be measured as a function of graph size rather than at a single point.
 *
 * Why that matters: the reflection providers replace native JVMCI metadata access, and a bigger
 * graph makes proportionally more metadata calls. So the overhead should SCALE with kernel
 * complexity, while the one-off JVM start-up difference does not. One kernel size cannot separate
 * those two; four sizes give a slope, and a slope predicts the cost for an arbitrary application.
 *
 * Three states are reported per tier:
 *   cold       first kernel of the run in a fresh JVM  - includes runtime start-up
 *   warm-new   every later kernel, each a NEW holder class, warm JIT
 *              -> this is the per-kernel compile cost, and the number that accumulates
 *   warm-same  re-execution of an already-compiled graph - everything cached, should be ~free
 *
 * Output is TSV on stdout (phase, tier, index, nanos); everything else goes to stderr.
 */
public final class CompilePerf {

    private static final int ELEMENTWISE_SIZE = 1024;
    private static final int DFT_SIZE = 256;    // O(n^2) kernel: keep device time small next to compile time
    private static final int NBODY_BODIES = 256; // ditto

    interface Kernel {
        /** Build the task graph for this kernel; the harness times its first execution. */
        void submit(TaskGraph tg, Workspace w);
    }

    /** Buffers for every tier, allocated once per graph so allocation never lands inside a timer. */
    static final class Workspace {
        final FloatArray a = new FloatArray(ELEMENTWISE_SIZE);
        final FloatArray b = new FloatArray(ELEMENTWISE_SIZE);
        final FloatArray c = new FloatArray(ELEMENTWISE_SIZE);
        final FloatArray inreal = new FloatArray(DFT_SIZE);
        final FloatArray inimag = new FloatArray(DFT_SIZE);
        final FloatArray outreal = new FloatArray(DFT_SIZE);
        final FloatArray outimag = new FloatArray(DFT_SIZE);
        final FloatArray pos = new FloatArray(NBODY_BODIES * 4);
        final FloatArray vel = new FloatArray(NBODY_BODIES * 4);

        Workspace() {
            a.init(2.0f);
            b.init(3.0f);
            c.init(0.0f);
            inreal.init(1.0f);
            inimag.init(0.5f);
            outreal.init(0.0f);
            outimag.init(0.0f);
            pos.init(1.0f);
            vel.init(0.0f);
        }
    }

    /** A kernel plus the tier it belongs to; the tier is what the results are grouped by. */
    record Entry(String tier, Kernel kernel) {
    }

    private static List<Entry> kernels() {
        List<Entry> ks = new ArrayList<>();
        // S: elementwise
        ks.add(new Entry("S", (tg, w) -> tg.task("t", Kernels.S0::run, w.a, w.b, w.c)));
        ks.add(new Entry("S", (tg, w) -> tg.task("t", Kernels.S1::run, w.a, w.b, w.c)));
        ks.add(new Entry("S", (tg, w) -> tg.task("t", Kernels.S2::run, w.a, w.b, w.c)));
        ks.add(new Entry("S", (tg, w) -> tg.task("t", Kernels.S3::run, w.a, w.b, w.c)));
        ks.add(new Entry("S", (tg, w) -> tg.task("t", Kernels.S4::run, w.a, w.b, w.c)));
        // M: math-heavy elementwise
        ks.add(new Entry("M", (tg, w) -> tg.task("t", Kernels.M0::run, w.a, w.b, w.c)));
        ks.add(new Entry("M", (tg, w) -> tg.task("t", Kernels.M1::run, w.a, w.b, w.c)));
        ks.add(new Entry("M", (tg, w) -> tg.task("t", Kernels.M2::run, w.a, w.b, w.c)));
        ks.add(new Entry("M", (tg, w) -> tg.task("t", Kernels.M3::run, w.a, w.b, w.c)));
        ks.add(new Entry("M", (tg, w) -> tg.task("t", Kernels.M4::run, w.a, w.b, w.c)));
        // L: DFT
        ks.add(new Entry("L", (tg, w) -> tg.task("t", Kernels.L0::run, w.inreal, w.inimag, w.outreal, w.outimag)));
        ks.add(new Entry("L", (tg, w) -> tg.task("t", Kernels.L1::run, w.inreal, w.inimag, w.outreal, w.outimag)));
        ks.add(new Entry("L", (tg, w) -> tg.task("t", Kernels.L2::run, w.inreal, w.inimag, w.outreal, w.outimag)));
        ks.add(new Entry("L", (tg, w) -> tg.task("t", Kernels.L3::run, w.inreal, w.inimag, w.outreal, w.outimag)));
        ks.add(new Entry("L", (tg, w) -> tg.task("t", Kernels.L4::run, w.inreal, w.inimag, w.outreal, w.outimag)));
        // XL: n-body
        ks.add(new Entry("XL", (tg, w) -> tg.task("t", Kernels.X0::run, NBODY_BODIES, w.pos, w.vel)));
        ks.add(new Entry("XL", (tg, w) -> tg.task("t", Kernels.X1::run, NBODY_BODIES, w.pos, w.vel)));
        ks.add(new Entry("XL", (tg, w) -> tg.task("t", Kernels.X2::run, NBODY_BODIES, w.pos, w.vel)));
        ks.add(new Entry("XL", (tg, w) -> tg.task("t", Kernels.X3::run, NBODY_BODIES, w.pos, w.vel)));
        ks.add(new Entry("XL", (tg, w) -> tg.task("t", Kernels.X4::run, NBODY_BODIES, w.pos, w.vel)));
        return ks;
    }

    private static long timeFirstExecution(Entry e, int graphId) throws Exception {
        Workspace w = new Workspace();
        // A distinct graph name per kernel: TornadoVM keys caches by task-graph identity, and
        // reusing one name would measure a graph update rather than a fresh compile.
        TaskGraph tg = new TaskGraph("probe" + graphId)
                .transferToDevice(DataTransferMode.EVERY_EXECUTION, w.a, w.b, w.inreal, w.inimag, w.pos, w.vel);
        e.kernel().submit(tg, w);
        tg.transferToHost(DataTransferMode.EVERY_EXECUTION, w.c, w.outreal, w.outimag, w.pos, w.vel);

        ImmutableTaskGraph itg = tg.snapshot();
        long start = System.nanoTime();
        try (TornadoExecutionPlan plan = new TornadoExecutionPlan(itg)) {
            plan.execute();
            long elapsed = System.nanoTime() - start;
            // Touch a result so nothing can be optimised away and a silently non-executing plan
            // shows up as a wrong value rather than as a suspiciously fast one.
            if (w.c.get(0) == Float.NEGATIVE_INFINITY) {
                System.err.println("unreachable " + w.outreal.get(0));
            }
            return elapsed;
        }
    }

    private static long[] timeWarmSame(Entry e, int graphId, int reps) throws Exception {
        Workspace w = new Workspace();
        TaskGraph tg = new TaskGraph("warm" + graphId)
                .transferToDevice(DataTransferMode.EVERY_EXECUTION, w.a, w.b, w.inreal, w.inimag, w.pos, w.vel);
        e.kernel().submit(tg, w);
        tg.transferToHost(DataTransferMode.EVERY_EXECUTION, w.c, w.outreal, w.outimag, w.pos, w.vel);

        long[] out = new long[reps];
        try (TornadoExecutionPlan plan = new TornadoExecutionPlan(tg.snapshot())) {
            plan.execute();   // first execution compiles; deliberately not recorded
            for (int i = 0; i < reps; i++) {
                long start = System.nanoTime();
                plan.execute();
                out[i] = System.nanoTime() - start;
            }
        }
        return out;
    }

    public static void main(String[] args) throws Exception {
        int warmReps = args.length > 0 ? Integer.parseInt(args[0]) : 20;
        List<Entry> ks = kernels();

        System.err.println("[probe] java=" + Runtime.version() + " kernels=" + ks.size() + " warmReps=" + warmReps);
        System.out.println("phase\ttier\tindex\tnanos");

        for (int i = 0; i < ks.size(); i++) {
            Entry e = ks.get(i);
            long ns = timeFirstExecution(e, i);
            // Only the very first kernel of the run is "cold": it carries runtime start-up. Every
            // later one is a warm-JIT, new-class compile, which is the per-kernel cost.
            System.out.println((i == 0 ? "cold" : "warm-new") + "\t" + e.tier() + "\t" + i + "\t" + ns);
            System.out.flush();
        }

        // Steady state, measured on the cheapest tier: if a per-execution regression existed it
        // would be most visible where there is least device work to hide it.
        for (long ns : timeWarmSame(ks.get(0), 0, warmReps)) {
            System.out.println("warm-same\tS\t0\t" + ns);
        }
        System.out.flush();
        System.err.println("[probe] done");
    }
}
