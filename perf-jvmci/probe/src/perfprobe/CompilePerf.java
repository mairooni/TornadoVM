package perfprobe;

import java.util.ArrayList;
import java.util.List;

import uk.ac.manchester.tornado.api.ImmutableTaskGraph;
import uk.ac.manchester.tornado.api.TaskGraph;
import uk.ac.manchester.tornado.api.TornadoExecutionPlan;
import uk.ac.manchester.tornado.api.enums.DataTransferMode;
import uk.ac.manchester.tornado.api.types.arrays.FloatArray;

/**
 * Prices TornadoVM's kernel-compilation path, so the JVMCI-removal overhead can be measured.
 *
 * TornadoVM compiles a task on its first execution, so "time to first execute" of a tiny kernel is
 * dominated by sketcher + Graal compile + code generation + driver compile. Three states are
 * reported, because "cold vs warm" hides the one that matters for a large application:
 *
 *   cold       first kernel in a fresh JVM. Pays classfile read + parse, cache fill, cold JIT.
 *   warm-new   kernels 1..N-1, each a NEW holder class, so the per-class metadata cost is paid
 *              again but the JIT is now warm. This is what an application with many kernels pays.
 *   warm-same  re-executions of kernel 0. Everything cached; should be ~free. If it is not, a
 *              cache is missing or being discarded.
 *
 * Output is TSV on stdout (phase, index, nanos) for 03-aggregate.py; everything else goes to
 * stderr so the two never mix.
 */
public final class CompilePerf {

    private static final int SIZE = 1024;

    interface Kernel {
        void submit(TaskGraph tg, FloatArray a, FloatArray b, FloatArray c);
    }

    /**
     * One entry per holder class. The method reference is what TornadoVM has to resolve back to a
     * java.lang.reflect.Method -- via JVMCI bytecode reading on the baseline branches, via
     * SerializedLambda + core reflection on the removal branch.
     */
    private static List<Kernel> kernels() {
        List<Kernel> ks = new ArrayList<>();
        ks.add((tg, a, b, c) -> tg.task("t", Kernels.K0::run, a, b, c));
        ks.add((tg, a, b, c) -> tg.task("t", Kernels.K1::run, a, b, c));
        ks.add((tg, a, b, c) -> tg.task("t", Kernels.K2::run, a, b, c));
        ks.add((tg, a, b, c) -> tg.task("t", Kernels.K3::run, a, b, c));
        ks.add((tg, a, b, c) -> tg.task("t", Kernels.K4::run, a, b, c));
        ks.add((tg, a, b, c) -> tg.task("t", Kernels.K5::run, a, b, c));
        ks.add((tg, a, b, c) -> tg.task("t", Kernels.K6::run, a, b, c));
        ks.add((tg, a, b, c) -> tg.task("t", Kernels.K7::run, a, b, c));
        return ks;
    }

    /** Build a fresh graph for `kernel` and time one execution end to end. */
    private static long timeFirstExecution(Kernel kernel, int graphId) throws Exception {
        FloatArray a = new FloatArray(SIZE);
        FloatArray b = new FloatArray(SIZE);
        FloatArray c = new FloatArray(SIZE);
        a.init(2.0f);
        b.init(3.0f);

        // A distinct graph name per kernel: TornadoVM keys some caches by task-graph identity, and
        // reusing one name across kernels would measure a graph update rather than a fresh compile.
        TaskGraph tg = new TaskGraph("probe" + graphId).transferToDevice(DataTransferMode.EVERY_EXECUTION, a, b);
        kernel.submit(tg, a, b, c);
        tg.transferToHost(DataTransferMode.EVERY_EXECUTION, c);

        ImmutableTaskGraph itg = tg.snapshot();
        long start = System.nanoTime();
        try (TornadoExecutionPlan plan = new TornadoExecutionPlan(itg)) {
            plan.execute();
            long elapsed = System.nanoTime() - start;
            // Touch the result so nothing above can be optimised away, and so a silently
            // non-executing plan shows up as a wrong value rather than a fast time.
            if (c.get(0) == Float.NEGATIVE_INFINITY) {
                System.err.println("unreachable " + c.get(0));
            }
            return elapsed;
        }
    }

    /** Re-execute one already-compiled graph: the steady-state, everything-cached path. */
    private static long[] timeWarmSame(Kernel kernel, int graphId, int reps) throws Exception {
        FloatArray a = new FloatArray(SIZE);
        FloatArray b = new FloatArray(SIZE);
        FloatArray c = new FloatArray(SIZE);
        a.init(2.0f);
        b.init(3.0f);

        TaskGraph tg = new TaskGraph("warm" + graphId).transferToDevice(DataTransferMode.EVERY_EXECUTION, a, b);
        kernel.submit(tg, a, b, c);
        tg.transferToHost(DataTransferMode.EVERY_EXECUTION, c);

        long[] out = new long[reps];
        try (TornadoExecutionPlan plan = new TornadoExecutionPlan(tg.snapshot())) {
            plan.execute();   // first execution compiles; not recorded here
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
        List<Kernel> ks = kernels();

        System.err.println("[probe] java=" + Runtime.version() + " kernels=" + ks.size() + " warmReps=" + warmReps);
        System.out.println("phase\tindex\tnanos");

        // cold = kernel 0; warm-new = kernels 1..n-1 (new holder class each, warm JIT)
        for (int i = 0; i < ks.size(); i++) {
            long ns = timeFirstExecution(ks.get(i), i);
            System.out.println((i == 0 ? "cold" : "warm-new") + "\t" + i + "\t" + ns);
            System.out.flush();
        }

        for (long ns : timeWarmSame(ks.get(0), 0, warmReps)) {
            System.out.println("warm-same\t0\t" + ns);
        }
        System.out.flush();
        System.err.println("[probe] done");
    }
}
