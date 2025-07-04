package uk.ac.manchester.tornado.examples;

import uk.ac.manchester.tornado.api.GridScheduler;
import uk.ac.manchester.tornado.api.ImmutableTaskGraph;
import uk.ac.manchester.tornado.api.KernelContext;
import uk.ac.manchester.tornado.api.TaskGraph;
import uk.ac.manchester.tornado.api.TornadoExecutionPlan;
import uk.ac.manchester.tornado.api.WorkerGrid;
import uk.ac.manchester.tornado.api.WorkerGrid1D;
import uk.ac.manchester.tornado.api.WorkerGrid2D;
import uk.ac.manchester.tornado.api.annotations.Parallel;
import uk.ac.manchester.tornado.api.enums.DataTransferMode;
import uk.ac.manchester.tornado.api.exceptions.TornadoExecutionPlanException;
import uk.ac.manchester.tornado.api.math.TornadoMath;
import uk.ac.manchester.tornado.api.types.matrix.Matrix2DDouble;

import java.util.ArrayList;
import java.util.LongSummaryStatistics;
import java.util.Random;

public class CalculateTVMDouble {
    // CHECKSTYLE:OFF
    // -- benchmark configuration
    public static final int WARM_UP_ITERATIONS = 10;
    public static final int MAX_ITERATIONS = 10;
    //---------------------------
    private static final double ORIGINAL_DEGREES_PER_RAD = 180.0d / Math.PI; // Approx 57.29577951308232
    private static final double ORIGINAL_EARTH_RADIUS_IN_METERS = ORIGINAL_DEGREES_PER_RAD * 60 * 1852; // Approx 6367468.459039645

    public static void calculate_tvm_kernel_race_free(Matrix2DDouble originalCoordsDegrees, Matrix2DDouble distances, int N) {

        final double DIVISOR_FOR_FROM_COORD = ORIGINAL_DEGREES_PER_RAD * ORIGINAL_DEGREES_PER_RAD;
        final double DIVISOR_FOR_TO_COORD = ORIGINAL_DEGREES_PER_RAD;

        for (@Parallel int from = 0; from < N - 1; from++) {
            double finalLatFrom = originalCoordsDegrees.get(from, 0) / DIVISOR_FOR_FROM_COORD;
            double finalLonFrom = originalCoordsDegrees.get(from, 1) / DIVISOR_FOR_FROM_COORD;

            for (int to = from + 1; to < N; to++) {
                double finalLatTo = originalCoordsDegrees.get(to, 0) / DIVISOR_FOR_TO_COORD;
                double finalLonTo = originalCoordsDegrees.get(to, 1) / DIVISOR_FOR_TO_COORD;

                double sinLatFrom = TornadoMath.sin(finalLatFrom);
                double sinLatTo = TornadoMath.sin(finalLatTo);
                double cosLatFrom = TornadoMath.cos(finalLatFrom);
                double cosLatTo = TornadoMath.cos(finalLatTo);
                double cosDeltaLon = TornadoMath.cos(finalLonFrom - finalLonTo);

                double argAcos = sinLatFrom * sinLatTo + cosLatFrom * cosLatTo * cosDeltaLon;

                argAcos = TornadoMath.max(-1.0, TornadoMath.min(1.0, argAcos));

                double acosResult = TornadoMath.acos(argAcos);
                double distance = acosResult * ORIGINAL_EARTH_RADIUS_IN_METERS;
                distances.set(from, to, distance);
            }
        }
    }

    public static void fillLowerTriangle(Matrix2DDouble distances, int N) {
        for (@Parallel int from = 1; from < N; from++) {
            for (int to = 0; to < from; to++) {
                distances.set(from, to, distances.get(to, from));
            }
        }
    }

    private static void calculate_tvm(double[][] coords, double[][] distances) {
        final double DEGREES_PER_RAD = 180.0d / Math.PI;
        final double EARTH_RADIUS_IN_METERS = DEGREES_PER_RAD * 60 * 1852;

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
                        distances[from][to] = Math.acos(Math.sin(coords[from][0]) * Math.sin(coords[to][0]) + Math.cos(coords[from][0]) * Math.cos(coords[to][0]) * Math.cos(coords[from][1] - coords[to][1])) * EARTH_RADIUS_IN_METERS;
                    } else {
                        distances[from][to] = 0;
                    }
                } else {
                    distances[from][to] = distances[to][from]; // !race condition: RW
                }
            }
        }
    }

    public static void setDiagonalToZero(Matrix2DDouble distances, int N) {
        for (@Parallel int i = 0; i < N; i++) {
            distances.set(i, i, 0.0);
        }
    }

    private static double[][] generateRandomCoordinates(int numPoints) {
        double[][] coords = new double[numPoints][2];
        Random rand = new Random(42);
        for (int i = 0; i < numPoints; i++) {
            double lat = -90 + 180 * rand.nextDouble();  // -90 to 90
            double lon = -180 + 360 * rand.nextDouble(); // -180 to 180
            coords[i][0] = lat;
            coords[i][1] = lon;
        }
        return coords;
    }

    public static void main(String[] args) throws Throwable {
        ArrayList<Integer> sizes = new ArrayList<>();
        int firstSize = 256;
        int numberOfSizes;

        if (args.length == 0) {
            numberOfSizes = 7;
        } else {
            numberOfSizes = Integer.parseInt(args[0]);
        }

        for (int i = 0; i < numberOfSizes; i++) {
            sizes.add(firstSize);
            firstSize *= 2;
        }
        for (Integer N : sizes) {

            System.out.println(">>> Device memory: " + System.getProperty("tornado.device.memory"));
            System.out.println("=== INPUT SIZE === " + N);
            ArrayList<Long> tornadoTimers = new ArrayList<>();
            ArrayList<Long> originalTimers = new ArrayList<>();

            double[][] coords = generateRandomCoordinates(N);

            // data for original java implementation
            double[][] coordsCopyOG = new double[N][2];
            for (int i = 0; i < N; i++) {
                coordsCopyOG[i][0] = coords[i][0];
                coordsCopyOG[i][1] = coords[i][1];
            }
            double[][] distancesSeqOG = new double[N][N];

            Matrix2DDouble inputCoordsMatrix = new Matrix2DDouble(N, 2);
            for (int i = 0; i < N; i++) {
                inputCoordsMatrix.set(i, 0, coords[i][0]);
                inputCoordsMatrix.set(i, 1, coords[i][1]);
            }

            Matrix2DDouble distancesTornado = new Matrix2DDouble(N, N);

            TaskGraph taskGraph = new TaskGraph("s0")
                    .transferToDevice(DataTransferMode.FIRST_EXECUTION, inputCoordsMatrix)
                    .task("t0", CalculateTVMDouble::setDiagonalToZero, distancesTornado, N)
                    .task("t1", CalculateTVMDouble::calculate_tvm_kernel_race_free, inputCoordsMatrix, distancesTornado, N)
                    .task("t2", CalculateTVMDouble::fillLowerTriangle, distancesTornado, N)
                    .transferToHost(DataTransferMode.EVERY_EXECUTION, distancesTornado);

            ImmutableTaskGraph immutableTaskGraph = taskGraph.snapshot();
            TornadoExecutionPlan executionPlan = new TornadoExecutionPlan(immutableTaskGraph);

            // -- tornado
            for (int i = 0; i < WARM_UP_ITERATIONS; i++) {
                executionPlan.execute();
            }

            for (int i = 0; i < MAX_ITERATIONS; i++) {
                long start = System.nanoTime();
                executionPlan.execute();
                long end = System.nanoTime();
                tornadoTimers.add((end - start));
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

            LongSummaryStatistics statsOG = originalTimers.stream().mapToLong(Long::longValue).summaryStatistics();
            LongSummaryStatistics statsTornado = tornadoTimers.stream().mapToLong(Long::longValue).summaryStatistics();

            System.out.println("OG    : " + statsOG.getAverage());
            System.out.println("TornadoVM: " + statsTornado.getAverage());
            System.out.println("SPEEDUP OVER OG: " + (statsOG.getAverage() / statsTornado.getAverage()));

            executionPlan.freeDeviceMemory();
        }
    }

    // CHECKSTYLE:ON
}
