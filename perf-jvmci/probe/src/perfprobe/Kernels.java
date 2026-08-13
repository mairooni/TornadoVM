package perfprobe;

import uk.ac.manchester.tornado.api.annotations.Parallel;
import uk.ac.manchester.tornado.api.math.TornadoMath;
import uk.ac.manchester.tornado.api.types.arrays.FloatArray;

/**
 * Kernel bodies for the compile-cost benchmark, in four complexity tiers.
 *
 * WHY TIERS. The cost we are pricing is metadata access during compilation -- resolving methods,
 * fields, types and constant-pool entries. A bigger graph makes more of those calls, so the
 * reflection overhead should grow with kernel complexity while the one-off JVM start-up saving does
 * not. Measuring one kernel size gives a single point that cannot distinguish those; measuring four
 * gives a slope, which is what lets anyone predict the cost for their own application.
 *
 * WHY ONE HOLDER CLASS PER KERNEL. The reflection path caches the parsed classfile per DECLARING
 * class (ReflectionUniverse.classfileBytes). Kernels sharing a holder would pay the read + parse
 * once and hide exactly the per-kernel cost being measured.
 *
 * The bodies are taken from tornado-benchmarks/ComputeKernels rather than invented, so the graph
 * shapes are the ones real TornadoVM applications actually compile:
 *   S  elementwise             - a few nodes, the floor
 *   M  reduction-ish with math - loop plus TornadoMath calls
 *   L  DFT                     - nested loop, trig, many locals
 *   XL n-body                  - triple-nested, array allocation, heaviest
 */
public final class Kernels {

    private Kernels() {
    }

    // ---------------------------------------------------------------- S: elementwise
    public static final class S0 {
        public static void run(FloatArray a, FloatArray b, FloatArray c) {
            for (@Parallel int i = 0; i < c.getSize(); i++) {
                c.set(i, a.get(i) + b.get(i));
            }
        }
    }

    public static final class S1 {
        public static void run(FloatArray a, FloatArray b, FloatArray c) {
            for (@Parallel int i = 0; i < c.getSize(); i++) {
                c.set(i, a.get(i) * b.get(i) + 1.0f);
            }
        }
    }

    public static final class S2 {
        public static void run(FloatArray a, FloatArray b, FloatArray c) {
            for (@Parallel int i = 0; i < c.getSize(); i++) {
                float v = a.get(i) - b.get(i);
                c.set(i, v * v);
            }
        }
    }

    public static final class S3 {
        public static void run(FloatArray a, FloatArray b, FloatArray c) {
            for (@Parallel int i = 0; i < c.getSize(); i++) {
                c.set(i, a.get(i) / (b.get(i) + 2.0f));
            }
        }
    }

    public static final class S4 {
        public static void run(FloatArray a, FloatArray b, FloatArray c) {
            for (@Parallel int i = 0; i < c.getSize(); i++) {
                c.set(i, (a.get(i) + b.get(i)) * 0.5f);
            }
        }
    }

    // ---------------------------------------------------------------- M: math-heavy elementwise
    public static final class M0 {
        public static void run(FloatArray a, FloatArray b, FloatArray c) {
            for (@Parallel int i = 0; i < c.getSize(); i++) {
                float v = TornadoMath.sqrt(a.get(i) * a.get(i) + b.get(i) * b.get(i));
                c.set(i, v + TornadoMath.exp(b.get(i) * 0.01f));
            }
        }
    }

    public static final class M1 {
        public static void run(FloatArray a, FloatArray b, FloatArray c) {
            for (@Parallel int i = 0; i < c.getSize(); i++) {
                float v = TornadoMath.sin(a.get(i)) * TornadoMath.cos(b.get(i));
                c.set(i, v * v + TornadoMath.abs(a.get(i)));
            }
        }
    }

    public static final class M2 {
        public static void run(FloatArray a, FloatArray b, FloatArray c) {
            for (@Parallel int i = 0; i < c.getSize(); i++) {
                float acc = 0.0f;
                for (int k = 0; k < 8; k++) {
                    acc += TornadoMath.sqrt(a.get(i) + k) * b.get(i);
                }
                c.set(i, acc);
            }
        }
    }

    public static final class M3 {
        public static void run(FloatArray a, FloatArray b, FloatArray c) {
            for (@Parallel int i = 0; i < c.getSize(); i++) {
                float acc = b.get(i);
                for (int k = 0; k < 8; k++) {
                    acc = acc * 0.5f + TornadoMath.exp(a.get(i) * 0.001f * k);
                }
                c.set(i, acc);
            }
        }
    }

    public static final class M4 {
        public static void run(FloatArray a, FloatArray b, FloatArray c) {
            for (@Parallel int i = 0; i < c.getSize(); i++) {
                float acc = 0.0f;
                for (int k = 0; k < 8; k++) {
                    acc += TornadoMath.sin(a.get(i) * k) + TornadoMath.cos(b.get(i) * k);
                }
                c.set(i, acc);
            }
        }
    }

    // ---------------------------------------------------------------- L: DFT (ComputeKernels.computeDFT)
    public static final class L0 {
        public static void run(FloatArray inreal, FloatArray inimag, FloatArray outreal, FloatArray outimag) {
            int n = inreal.getSize();
            for (@Parallel int k = 0; k < n; k++) {
                float sumreal = 0;
                float sumimag = 0;
                for (int t = 0; t < n; t++) {
                    float angle = ((2 * TornadoMath.floatPI() * t * k) / n);
                    sumreal += (inreal.get(t) * TornadoMath.cos(angle)) + (inimag.get(t) * TornadoMath.sin(angle));
                    sumimag += -(inreal.get(t) * TornadoMath.sin(angle)) + (inimag.get(t) * TornadoMath.cos(angle));
                }
                outreal.set(k, sumreal);
                outimag.set(k, sumimag);
            }
        }
    }

    public static final class L1 {
        public static void run(FloatArray inreal, FloatArray inimag, FloatArray outreal, FloatArray outimag) {
            int n = inreal.getSize();
            for (@Parallel int k = 0; k < n; k++) {
                float sumreal = 0;
                float sumimag = 0;
                for (int t = 0; t < n; t++) {
                    float angle = ((2 * TornadoMath.floatPI() * t * k) / n) + 0.5f;
                    sumreal += (inreal.get(t) * TornadoMath.cos(angle)) - (inimag.get(t) * TornadoMath.sin(angle));
                    sumimag += (inreal.get(t) * TornadoMath.sin(angle)) + (inimag.get(t) * TornadoMath.cos(angle));
                }
                outreal.set(k, sumreal);
                outimag.set(k, sumimag);
            }
        }
    }

    public static final class L2 {
        public static void run(FloatArray inreal, FloatArray inimag, FloatArray outreal, FloatArray outimag) {
            int n = inreal.getSize();
            for (@Parallel int k = 0; k < n; k++) {
                float sumreal = 0;
                float sumimag = 0;
                for (int t = 0; t < n; t++) {
                    float angle = ((2 * TornadoMath.floatPI() * t * k) / n);
                    float w = TornadoMath.sqrt(1.0f + t);
                    sumreal += w * (inreal.get(t) * TornadoMath.cos(angle));
                    sumimag += w * (inimag.get(t) * TornadoMath.sin(angle));
                }
                outreal.set(k, sumreal);
                outimag.set(k, sumimag);
            }
        }
    }

    public static final class L3 {
        public static void run(FloatArray inreal, FloatArray inimag, FloatArray outreal, FloatArray outimag) {
            int n = inreal.getSize();
            for (@Parallel int k = 0; k < n; k++) {
                float sumreal = 0;
                float sumimag = 0;
                for (int t = 0; t < n; t++) {
                    float angle = ((2 * TornadoMath.floatPI() * t * k) / n) * 0.5f;
                    sumreal += (inreal.get(t) * TornadoMath.cos(angle)) + TornadoMath.exp(-0.001f * t);
                    sumimag += (inimag.get(t) * TornadoMath.sin(angle));
                }
                outreal.set(k, sumreal);
                outimag.set(k, sumimag);
            }
        }
    }

    public static final class L4 {
        public static void run(FloatArray inreal, FloatArray inimag, FloatArray outreal, FloatArray outimag) {
            int n = inreal.getSize();
            for (@Parallel int k = 0; k < n; k++) {
                float sumreal = 0;
                float sumimag = 0;
                for (int t = 0; t < n; t++) {
                    float angle = ((2 * TornadoMath.floatPI() * t * k) / n);
                    sumreal += (inreal.get(t) * TornadoMath.cos(angle)) * (1.0f + TornadoMath.abs(inimag.get(t)));
                    sumimag += (inimag.get(t) * TornadoMath.sin(angle)) * (1.0f + TornadoMath.abs(inreal.get(t)));
                }
                outreal.set(k, sumreal);
                outimag.set(k, sumimag);
            }
        }
    }

    // ---------------------------------------------------------------- XL: n-body (ComputeKernels.nBody)
    public static final class X0 {
        public static void run(int numBodies, FloatArray refPos, FloatArray refVel) {
            final float delT = 0.005f;
            final float espSqr = 500.0f;
            for (@Parallel int i = 0; i < numBodies; i++) {
                int body = 4 * i;
                float[] acc = new float[] { 0.0f, 0.0f, 0.0f };
                for (int j = 0; j < numBodies; j++) {
                    float[] r = new float[3];
                    int index = 4 * j;
                    float distSqr = 0.0f;
                    for (int k = 0; k < 3; k++) {
                        r[k] = refPos.get(index + k) - refPos.get(body + k);
                        distSqr += r[k] * r[k];
                    }
                    float invDist = 1.0f / TornadoMath.sqrt(distSqr + espSqr);
                    float invDistCube = invDist * invDist * invDist;
                    float s = refPos.get(index + 3) * invDistCube;
                    for (int k = 0; k < 3; k++) {
                        acc[k] += s * r[k];
                    }
                }
                for (int k = 0; k < 3; k++) {
                    refPos.set(body + k, refPos.get(body + k) + refPos.get(body + k) * delT + 0.5f * acc[k] * delT * delT);
                    refVel.set(body + k, refPos.get(body + k) + acc[k] * delT);
                }
            }
        }
    }

    public static final class X1 {
        public static void run(int numBodies, FloatArray refPos, FloatArray refVel) {
            final float delT = 0.004f;
            final float espSqr = 400.0f;
            for (@Parallel int i = 0; i < numBodies; i++) {
                int body = 4 * i;
                float[] acc = new float[] { 0.0f, 0.0f, 0.0f };
                for (int j = 0; j < numBodies; j++) {
                    float[] r = new float[3];
                    int index = 4 * j;
                    float distSqr = 0.0f;
                    for (int k = 0; k < 3; k++) {
                        r[k] = refPos.get(index + k) - refPos.get(body + k);
                        distSqr += r[k] * r[k];
                    }
                    float invDist = 1.0f / TornadoMath.sqrt(distSqr + espSqr);
                    float s = refPos.get(index + 3) * invDist * invDist * invDist;
                    for (int k = 0; k < 3; k++) {
                        acc[k] += s * r[k] * 0.5f;
                    }
                }
                for (int k = 0; k < 3; k++) {
                    refPos.set(body + k, refPos.get(body + k) + acc[k] * delT * delT);
                    refVel.set(body + k, refVel.get(body + k) + acc[k] * delT);
                }
            }
        }
    }

    public static final class X2 {
        public static void run(int numBodies, FloatArray refPos, FloatArray refVel) {
            final float delT = 0.006f;
            final float espSqr = 600.0f;
            for (@Parallel int i = 0; i < numBodies; i++) {
                int body = 4 * i;
                float[] acc = new float[] { 0.0f, 0.0f, 0.0f };
                for (int j = 0; j < numBodies; j++) {
                    float[] r = new float[3];
                    int index = 4 * j;
                    float distSqr = 0.0f;
                    for (int k = 0; k < 3; k++) {
                        r[k] = refPos.get(index + k) - refPos.get(body + k);
                        distSqr += r[k] * r[k];
                    }
                    float invDist = 1.0f / TornadoMath.sqrt(distSqr + espSqr);
                    float invDistCube = invDist * invDist * invDist;
                    float s = refPos.get(index + 3) * invDistCube;
                    for (int k = 0; k < 3; k++) {
                        acc[k] += s * r[k] + 0.001f;
                    }
                }
                for (int k = 0; k < 3; k++) {
                    refPos.set(body + k, refPos.get(body + k) + refPos.get(body + k) * delT);
                    refVel.set(body + k, refPos.get(body + k) * 0.5f + acc[k] * delT);
                }
            }
        }
    }

    public static final class X3 {
        public static void run(int numBodies, FloatArray refPos, FloatArray refVel) {
            final float delT = 0.003f;
            final float espSqr = 300.0f;
            for (@Parallel int i = 0; i < numBodies; i++) {
                int body = 4 * i;
                float[] acc = new float[] { 0.0f, 0.0f, 0.0f };
                for (int j = 0; j < numBodies; j++) {
                    float[] r = new float[3];
                    int index = 4 * j;
                    float distSqr = 0.0f;
                    for (int k = 0; k < 3; k++) {
                        r[k] = refPos.get(index + k) - refPos.get(body + k);
                        distSqr += r[k] * r[k] * 1.0001f;
                    }
                    float invDist = 1.0f / TornadoMath.sqrt(distSqr + espSqr);
                    float s = refPos.get(index + 3) * invDist * invDist * invDist;
                    for (int k = 0; k < 3; k++) {
                        acc[k] += s * r[k];
                    }
                }
                for (int k = 0; k < 3; k++) {
                    refPos.set(body + k, refPos.get(body + k) + acc[k] * delT * delT * 0.5f);
                    refVel.set(body + k, refVel.get(body + k) + acc[k] * delT);
                }
            }
        }
    }

    public static final class X4 {
        public static void run(int numBodies, FloatArray refPos, FloatArray refVel) {
            final float delT = 0.007f;
            final float espSqr = 700.0f;
            for (@Parallel int i = 0; i < numBodies; i++) {
                int body = 4 * i;
                float[] acc = new float[] { 0.0f, 0.0f, 0.0f };
                for (int j = 0; j < numBodies; j++) {
                    float[] r = new float[3];
                    int index = 4 * j;
                    float distSqr = 0.0f;
                    for (int k = 0; k < 3; k++) {
                        r[k] = refPos.get(index + k) - refPos.get(body + k);
                        distSqr += r[k] * r[k];
                    }
                    float invDist = 1.0f / TornadoMath.sqrt(distSqr + espSqr + 1.0f);
                    float invDistCube = invDist * invDist * invDist;
                    float s = refPos.get(index + 3) * invDistCube;
                    for (int k = 0; k < 3; k++) {
                        acc[k] += s * r[k];
                    }
                }
                for (int k = 0; k < 3; k++) {
                    refPos.set(body + k, refPos.get(body + k) * 1.0001f + acc[k] * delT);
                    refVel.set(body + k, refVel.get(body + k) + acc[k] * delT * 0.5f);
                }
            }
        }
    }
}
