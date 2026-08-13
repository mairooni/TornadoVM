package perfprobe;

import uk.ac.manchester.tornado.api.annotations.Parallel;
import uk.ac.manchester.tornado.api.types.arrays.FloatArray;

/**
 * Kernel bodies for the compile-cost probe, spread over several holder classes.
 *
 * Several classes, not one: the reflection metadata path caches the parsed classfile per
 * DECLARING CLASS (ReflectionUniverse.classfileBytes), so kernels sharing a holder would pay the
 * classfile read + parse once and hide exactly the cost we are trying to price. One kernel per
 * holder makes every "warm-new" kernel a genuine new-class miss.
 *
 * The bodies are deliberately trivial and the arrays small: we are timing COMPILATION, so device
 * execution and transfers must stay negligible next to it.
 */
public final class Kernels {

    private Kernels() {
    }

    public static final class K0 {
        public static void run(FloatArray a, FloatArray b, FloatArray c) {
            for (@Parallel int i = 0; i < c.getSize(); i++) {
                c.set(i, a.get(i) + b.get(i));
            }
        }
    }

    public static final class K1 {
        public static void run(FloatArray a, FloatArray b, FloatArray c) {
            for (@Parallel int i = 0; i < c.getSize(); i++) {
                c.set(i, a.get(i) * b.get(i) + 1.0f);
            }
        }
    }

    public static final class K2 {
        public static void run(FloatArray a, FloatArray b, FloatArray c) {
            for (@Parallel int i = 0; i < c.getSize(); i++) {
                float v = a.get(i) - b.get(i);
                c.set(i, v * v);
            }
        }
    }

    public static final class K3 {
        public static void run(FloatArray a, FloatArray b, FloatArray c) {
            for (@Parallel int i = 0; i < c.getSize(); i++) {
                c.set(i, a.get(i) / (b.get(i) + 2.0f));
            }
        }
    }

    public static final class K4 {
        public static void run(FloatArray a, FloatArray b, FloatArray c) {
            for (@Parallel int i = 0; i < c.getSize(); i++) {
                c.set(i, a.get(i) + b.get(i) * 3.0f);
            }
        }
    }

    public static final class K5 {
        public static void run(FloatArray a, FloatArray b, FloatArray c) {
            for (@Parallel int i = 0; i < c.getSize(); i++) {
                float v = a.get(i) + 1.0f;
                c.set(i, v * b.get(i));
            }
        }
    }

    public static final class K6 {
        public static void run(FloatArray a, FloatArray b, FloatArray c) {
            for (@Parallel int i = 0; i < c.getSize(); i++) {
                c.set(i, (a.get(i) + b.get(i)) * 0.5f);
            }
        }
    }

    public static final class K7 {
        public static void run(FloatArray a, FloatArray b, FloatArray c) {
            for (@Parallel int i = 0; i < c.getSize(); i++) {
                float v = b.get(i) - 1.0f;
                c.set(i, a.get(i) + v);
            }
        }
    }
}
