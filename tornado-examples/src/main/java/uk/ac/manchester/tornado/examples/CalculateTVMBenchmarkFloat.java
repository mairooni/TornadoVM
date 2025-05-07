package uk.ac.manchester.tornado.examples;

import uk.ac.manchester.tornado.api.TaskGraph;
import uk.ac.manchester.tornado.api.TornadoExecutionPlan;
import uk.ac.manchester.tornado.api.annotations.Parallel;
import uk.ac.manchester.tornado.api.enums.DataTransferMode;
import uk.ac.manchester.tornado.api.math.TornadoMath;
import uk.ac.manchester.tornado.api.types.matrix.Matrix2DFloat;

import java.util.ArrayList;
import java.util.LongSummaryStatistics;
import java.util.Random;

public class CalculateTVMBenchmarkFloat {
    // CHECKSTYLE:OFF
    public static final int WARM_UP_ITERATIONS = 10;
    public static final int MAX_ITERATIONS = 10;

    // tornadovm compatible implementation
    public static void normalizeCoords(Matrix2DFloat coords, int N) {
        for (@Parallel int i = 0; i < N; i++) {
            coords.set(i, 0, coords.get(i, 0) / (float) (180 / Math.PI));
            coords.set(i, 1, coords.get(i, 1) / (float) (180 / Math.PI));
        }
    }

    public static void computeUpperTriangle(Matrix2DFloat coords, Matrix2DFloat distances, int N) {
        final float EARTH_RADIUS_IN_METERS = (float) ((180 / Math.PI) * 60 * 1852);

        for (@Parallel int from = 0; from < N; from++) {
            for (int to = from; to < N; to++) {
                float distance;
                if (from == to) {
                    distance = 0;
                } else {
                    float lat1 = coords.get(from, 0);
                    float lon1 = coords.get(from, 1);
                    float lat2 = coords.get(to, 0);
                    float lon2 = coords.get(to, 1);

                    distance = TornadoMath.acos(
                            TornadoMath.sin(lat1) * TornadoMath.sin(lat2) +
                                    TornadoMath.cos(lat1) * TornadoMath.cos(lat2) * TornadoMath.cos(lon1 - lon2)
                    ) * EARTH_RADIUS_IN_METERS;
                }
                distances.set(from, to, distance);
            }
        }
    }

    public static void mirrorLowerTriangle(Matrix2DFloat distances, int N) {
        for (@Parallel int from = 0; from < N; from++) {
            for (int to = 0; to < from; to++) {
                distances.set(from, to, distances.get(to, from));
            }
        }
    }

    // original implementation
    private static void calculate_tvm(float[][] coords, float[][] distances) {
        final float DEGREES_PER_RAD = (float) (180 / Math.PI);
        final float EARTH_RADIUS_IN_METERS = DEGREES_PER_RAD * 60 * 1852;

        for (int from = 0; from < coords.length; from++) {
            coords[from][0] /= DEGREES_PER_RAD;
            coords[from][1] /= DEGREES_PER_RAD;
            for (int to = 0; to < coords.length; to++) {
                if (from <= to) {
                    if (from == 0) { // !race condition:
                        coords[to][0] /= DEGREES_PER_RAD;
                        coords[to][1] /= DEGREES_PER_RAD;
                    }
                    if (from < to) {
                        distances[from][to] = (float) (Math.acos(Math.sin(coords[from][0]) * Math.sin(coords[to][0]) + Math.cos(coords[from][0]) * Math.cos(coords[to][0]) * Math.cos(coords[from][1] - coords[to][1])) * EARTH_RADIUS_IN_METERS);
                    } else {
                        distances[from][to] = 0;
                    }
                } else {
                    distances[from][to] = distances[to][from]; // !race condition: RW
                }
            }
        }
    }

    // data generator
    private static float[][] generateRandomCoordinates(int numPoints) {
        float[][] coords = new float[numPoints][2];
        Random rand = new Random(42);
        for (int i = 0; i < numPoints; i++) {
            float lat = -90 + 180 * rand.nextFloat();  // -90 to 90
            float lon = -180 + 360 * rand.nextFloat(); // -180 to 180
            coords[i][0] = lat;
            coords[i][1] = lon;
        }
        return coords;
    }

    public static void main(String args[]) throws Throwable {
        ArrayList<Integer> sizes = new ArrayList<>();
        int firstSize = 64;
        int numberOfSizes;

        if (args.length == 0) {
            numberOfSizes = 24;
        } else {
            numberOfSizes = Integer.parseInt(args[0]);
        }

        for (int i = 0; i < numberOfSizes; i++) {
            sizes.add(firstSize);
            firstSize *= 2;
        }

        for (Integer N : sizes) {
            System.out.println("=== INPUT SIZE === " + N);
            ArrayList<Long> seqTimers = new ArrayList<>();
            ArrayList<Long> tornadoTimers = new ArrayList<>();
            ArrayList<Long> originalTimers = new ArrayList<>();

            // data for tornadovm version
            float[][] coords = generateRandomCoordinates(N);
            Matrix2DFloat coordsM = new Matrix2DFloat(coords);
            float[][] distances = new float[N][N];
            Matrix2DFloat distancesM = new Matrix2DFloat(distances);
            // data for tornadovm java version
            float[][] coordsCopy = new float[N][2];
            for (int i = 0; i < N; i++) {
                coordsCopy[i][0] = coords[i][0];
                coordsCopy[i][1] = coords[i][1];
            }
            Matrix2DFloat copySeqM = new Matrix2DFloat(coordsCopy);
            float[][] distancesSeq = new float[N][N];
            Matrix2DFloat distancesSeqM = new Matrix2DFloat(distancesSeq);
            // data for original java implementation
            float[][] coordsCopyOG = new float[N][2];
            for (int i = 0; i < N; i++) {
                coordsCopyOG[i][0] = coords[i][0];
                coordsCopyOG[i][1] = coords[i][1];
            }
            float[][] distancesSeqOG = new float[N][N];

            // TornadoVM task graph declaration
            TaskGraph taskGraph = new TaskGraph("s0")
                    .transferToDevice(DataTransferMode.FIRST_EXECUTION, coordsM, distancesM)
                    .task("t0", CalculateTVMBenchmarkFloat::normalizeCoords, coordsM, N)
                    .task("t1", CalculateTVMBenchmarkFloat::computeUpperTriangle, coordsM, distancesM, N)
                    .task("t2", CalculateTVMBenchmarkFloat::mirrorLowerTriangle, distancesM, N)
                    .transferToHost(DataTransferMode.EVERY_EXECUTION, distancesM);

            TornadoExecutionPlan executor = new TornadoExecutionPlan(taskGraph.snapshot());

            // -- tornado
            for (int i = 0; i < WARM_UP_ITERATIONS; i++) {
                executor.execute();
            }

            for (int i = 0; i < MAX_ITERATIONS; i++) {
                long start = System.nanoTime();
                executor.execute();
                long end = System.nanoTime();
                tornadoTimers.add((end - start));
            }

            // -- sequential
            for (int i = 0; i < WARM_UP_ITERATIONS; i++) {
                normalizeCoords(copySeqM, N);
                computeUpperTriangle(copySeqM, distancesSeqM, N);
                mirrorLowerTriangle(distancesSeqM, N);
            }

            for (int i = 0; i < MAX_ITERATIONS; i++) {
                long start = System.nanoTime();
                normalizeCoords(copySeqM, N);
                computeUpperTriangle(copySeqM, distancesSeqM, N);
                mirrorLowerTriangle(distancesSeqM, N);
                long end = System.nanoTime();
                seqTimers.add((end - start));
            }

            // -- original
            for (int i = 0; i < WARM_UP_ITERATIONS; i++) {
                calculate_tvm(coordsCopyOG, distancesSeqOG);
            }

            for (int i = 0; i < MAX_ITERATIONS; i++) {
                long start = System.nanoTime();
                calculate_tvm(coordsCopyOG, distancesSeqOG);
                long end = System.nanoTime();
                originalTimers.add((end - start));
            }

            // Compute Medians
            LongSummaryStatistics statsSeq = seqTimers.stream().mapToLong(Long::longValue).summaryStatistics();
            LongSummaryStatistics statsOG = originalTimers.stream().mapToLong(Long::longValue).summaryStatistics();
            LongSummaryStatistics statsTornado = tornadoTimers.stream().mapToLong(Long::longValue).summaryStatistics();

            System.out.println("SEQ    : " + statsSeq.getAverage());
            System.out.println("OG    : " + statsOG.getAverage());
            System.out.println("TornadoVM: " + statsTornado.getAverage());
            System.out.println("SPEEDUP OVER SEQ: " + (statsSeq.getAverage() / statsTornado.getAverage()));
            System.out.println("SPEEDUP OVER OG: " + (statsOG.getAverage() / statsTornado.getAverage()));
            System.out.println("SPEEDUP SEQ | SPEEDUP OG");
            System.out.println((statsSeq.getAverage() / statsTornado.getAverage()) + " " + (statsOG.getAverage() / statsTornado.getAverage()));

            executor.freeDeviceMemory();
        }
    }

    // CHECKSTYLE:ON
}
