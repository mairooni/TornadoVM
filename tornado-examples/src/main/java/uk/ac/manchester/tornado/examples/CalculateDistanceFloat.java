package uk.ac.manchester.tornado.examples;

import uk.ac.manchester.tornado.api.ImmutableTaskGraph;
import uk.ac.manchester.tornado.api.TaskGraph;
import uk.ac.manchester.tornado.api.TornadoExecutionPlan;
import uk.ac.manchester.tornado.api.annotations.Parallel;
import uk.ac.manchester.tornado.api.enums.DataTransferMode;
import uk.ac.manchester.tornado.api.math.TornadoMath;

import java.util.ArrayList;
import java.util.LongSummaryStatistics;
import java.util.Random;

public class CalculateDistanceFloat {
    // CHECKSTYLE:OFF

    private static float EARTH_RADIUS = 6378;
    public static final int WARM_UP_ITERATIONS = 10;
    public static final int MAX_ITERATIONS = 10;

    public static void calculateDistance(float[] lat1, float[] lon1, float[] lat2, float[] lon2, float[] distance) {
        for (@Parallel int i = 0; i < lat1.length; i++) {
            float lat1Rad = TornadoMath.toRadians(lat1[i]);
            float lon1Rad = TornadoMath.toRadians(lon1[i]);
            float lat2Rad = TornadoMath.toRadians(lat2[i]);
            float lon2Rad = TornadoMath.toRadians(lon2[i]);

            float x = (lon2Rad - lon1Rad) * TornadoMath.cos((lat1Rad + lat2Rad) / 2);
            float y = (lat2Rad - lat1Rad);
            distance[i] = TornadoMath.sqrt(x * x + y * y) * EARTH_RADIUS;
        }
    }

    public static void main(String[] args) throws Throwable {
        ArrayList<Integer> sizes = new ArrayList<>();
        int firstSize = 32768;
        int numberOfSizes;

        if (args.length == 0) {
            numberOfSizes = 12;
        } else {
            numberOfSizes = Integer.parseInt(args[0]);
        }

        for (int i = 0; i < numberOfSizes; i++) {
            sizes.add(firstSize);
            firstSize *= 2;
        }

        for (Integer N : sizes) {
            System.out.println("=== INPUT SIZE === " + N);
            float[] lat1 = new float[N];
            float[] lon1 = new float[N];
            float[] lat2 = new float[N];
            float[] lon2 = new float[N];
            float[] distance = new float[N];
            float[] distanceSeq = new float[N];

            Random random = new Random();

            for (int i = 0; i < N; i++) {
                // Random lat: -90 to +90
                lat1[i] = -90 + 180 * random.nextFloat();
                lat2[i] = -90 + 180 * random.nextFloat();

                // Random lon: -180 to +180
                lon1[i] = -180 + 360 * random.nextFloat();
                lon2[i] = -180 + 360 * random.nextFloat();
            }

            ArrayList<Long> seqTimers = new ArrayList<>();
            ArrayList<Long> tornadoTimers = new ArrayList<>();

            TaskGraph taskGraph = new TaskGraph("s0") //
                    .transferToDevice(DataTransferMode.FIRST_EXECUTION, lat1, lon1, lat2, lon2) //
                    .task("t0", CalculateDistanceFloat::calculateDistance, lat1, lon1, lat2, lon2, distance) //
                    .transferToHost(DataTransferMode.EVERY_EXECUTION, distance);

            ImmutableTaskGraph immutableTaskGraph = taskGraph.snapshot();
            TornadoExecutionPlan executor = new TornadoExecutionPlan(immutableTaskGraph);

            // -- sequential
            for (int i = 0; i < WARM_UP_ITERATIONS; i++) {
                calculateDistance(lat1, lon1, lat2, lon2, distanceSeq);
            }

            for (int i = 0; i < MAX_ITERATIONS; i++) {
                long start = System.nanoTime();
                calculateDistance(lat1, lon1, lat2, lon2, distanceSeq);
                long end = System.nanoTime();
                seqTimers.add((end - start));
            }

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

            // Compute Medians
            LongSummaryStatistics statsSeq = seqTimers.stream().mapToLong(Long::longValue).summaryStatistics();
            LongSummaryStatistics statsTornado = tornadoTimers.stream().mapToLong(Long::longValue).summaryStatistics();

            System.out.println("SEQ    : " + statsSeq.getAverage());
            System.out.println("TornadoVM: " + statsTornado.getAverage());
            System.out.println("SPEEDUP: " + (statsSeq.getAverage() / statsTornado.getAverage()));

            executor.freeDeviceMemory();
        }

    }

    // CHECKSTYLE:ON
}
