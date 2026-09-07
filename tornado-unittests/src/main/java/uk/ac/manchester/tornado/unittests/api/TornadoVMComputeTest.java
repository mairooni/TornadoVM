package uk.ac.manchester.tornado.unittests.api;

import uk.ac.manchester.tornado.api.AccessorParameters;
import uk.ac.manchester.tornado.api.TaskGraph;
import uk.ac.manchester.tornado.api.TornadoExecutionPlan;
import uk.ac.manchester.tornado.api.ImmutableTaskGraph;
import uk.ac.manchester.tornado.api.common.Access;
import uk.ac.manchester.tornado.api.enums.DataTransferMode;
import uk.ac.manchester.tornado.api.types.arrays.FloatArray;
import uk.ac.manchester.tornado.api.types.arrays.IntArray;
import uk.ac.manchester.tornado.api.annotations.Parallel;

import org.junit.Test;
import uk.ac.manchester.tornado.unittests.common.TornadoTestBase;

import static org.junit.Assert.assertTrue;

/**
 * Unit test for testing a multi-head attention mechanism on TornadoVM
 * with specific dimensions:
 * dim = 2048, headSize = 64, numHeads = 32, numKVHeads = 8
 * Run uk.ac.manchester.tornado.unittests.api.TornadoVMComputeTest#testGroupedQueryAttention
 */
public class TornadoVMComputeTest extends TornadoTestBase {

    public static class TornadoVMCompute {


        public static void processHeadsParallel(
                FloatArray q, FloatArray key_cache, FloatArray value_cache, FloatArray xb,
                int nHeads, int headSize, int kvDim, int kvMul, int seqLen,
                IntArray positionNlayer, FloatArray wrapAtt) {

            int pos = positionNlayer.get(0);
            int layer = positionNlayer.get(1);
            long loff = layer * seqLen * kvDim; // layer offset into KV cache

            // Parallelize computation across attention heads
            for (@Parallel int h = 0; h < nHeads; h++) {
                // Process each head in parallel
                processHeadTornado(q, key_cache, value_cache, xb, h, headSize, kvDim, kvMul, loff, pos, wrapAtt);
            }
        }

        private static void processHeadTornado(
                FloatArray allQ, FloatArray key_cache, FloatArray value_cache, FloatArray allXb,
                int h, int headSize, int kvDim, int kvMul, long loff, int pos, FloatArray wrapAtt) {

            // Base index for this head's attention weights
            int headOffset = h * (pos + 1);

            // STEP 1: Calculate attention scores for all timesteps
            for (int t = 0; t <= pos; t++) {
                int kvHeadIdx = h / kvMul;
                int keyOffset = (int) (loff + t * kvDim + kvHeadIdx * headSize);

                float score = 0.0f;
                for (int i = 0; i < headSize; i++) {
                    score += allQ.get(h * headSize + i) * key_cache.get(keyOffset + i);
                }
                score = score / (float) Math.sqrt(headSize);

                // Store in attention buffer
                wrapAtt.set(headOffset + t, score);
            }

            // STEP 2: Find max score for softmax stability
            float maxScore = wrapAtt.get(headOffset);
            for (int t = 1; t <= pos; t++) {
                float val = wrapAtt.get(headOffset + t);
                if (val > maxScore) {
                    maxScore = val;
                }
            }

            // STEP 3: Compute exponentials and sum
            float sum = 0.0f;
            for (int t = 0; t <= pos; t++) {
                int idx = headOffset + t;
                float expScore = (float) Math.exp(wrapAtt.get(idx) - maxScore);
                wrapAtt.set(idx, expScore);
                sum += expScore;
            }

            // STEP 4: Normalize
            float normFactor = (sum > 0.0f) ? (1.0f / sum) : (1.0f / (pos + 1));
            for (int t = 0; t <= pos; t++) {
                int idx = headOffset + t;
                wrapAtt.set(idx, wrapAtt.get(idx) * normFactor);
            }

            // STEP 5: Compute weighted sum of values for each dimension
            for (int i = 0; i < headSize; i++) {
                float weightedSum = 0.0f;
                for (int t = 0; t <= pos; t++) {
                    int kvHeadIdx = h / kvMul;
                    int valueOffset = (int) (loff + t * kvDim + kvHeadIdx * headSize);
                    weightedSum += wrapAtt.get(headOffset + t) * value_cache.get(valueOffset + i);
                }
                allXb.set(h * headSize + i, weightedSum);
            }
        }
    }

    /**
     * Sequential reference implementation of multi-head attention
     * Used to validate the parallel TornadoVM implementation
     */
    private static void processAttentionSequential(
            FloatArray query, FloatArray keyCache, FloatArray valueCache, FloatArray output,
            int nHeads, int headSize, int kvDim, int kvMul, int seqLen, int pos, int layer) {

        // Process each attention head
        for (int h = 0; h < nHeads; h++) {
            int layerOffset = layer * seqLen * kvDim;

            // Local scores for this head
            float[] scores = new float[pos + 1];

            // Calculate attention scores
            for (int t = 0; t <= pos; t++) {
                int kvHeadIndex = h / kvMul;
                int keyOffset = layerOffset + t * kvDim + kvHeadIndex * headSize;

                float score = 0.0f;
                for (int i = 0; i < headSize; i++) {
                    score += query.get(h * headSize + i) * keyCache.get(keyOffset + i);
                }

                // Scale by sqrt(headSize)
                score = score / (float)Math.sqrt(headSize);
                scores[t] = score;
            }

            // Apply softmax
            float maxScore = Float.NEGATIVE_INFINITY;
            for (int t = 0; t <= pos; t++) {
                if (scores[t] > maxScore) {
                    maxScore = scores[t];
                }
            }

            float sum = 0.0f;
            for (int t = 0; t <= pos; t++) {
                scores[t] = (float)Math.exp(scores[t] - maxScore);
                sum += scores[t];
            }

            for (int t = 0; t <= pos; t++) {
                scores[t] = scores[t] / sum;
            }

            // Compute weighted sum of values
            for (int i = 0; i < headSize; i++) {
                float weightedSum = 0.0f;

                for (int t = 0; t <= pos; t++) {
                    int kvHeadIndex = h / kvMul;
                    int valueOffset = layerOffset + t * kvDim + kvHeadIndex * headSize + i;
                    weightedSum += scores[t] * valueCache.get(valueOffset);
                }

                // Store in output
                output.set(h * headSize + i, weightedSum);
            }
        }
    }

    /**
     * Test to verify multi-head attention implementation using TornadoVM
     * with dimensions: dim = 2048, headSize = 64, numHeads = 32, numKVHeads = 8
     */
    @Test
    public void testGroupedQueryAttention() {
        // Model dimensions and parameters based on the specification
        final int dim = 2048;
        final int nHeads = 32;
        final int headSize = 64;
        final int numKVHeads = 8;
        final int kvMul = nHeads / numKVHeads; // Grouped Query Attention factor
        final int kvDim = headSize * numKVHeads;
        final int seqLen = 128;
        final int pos = 16;       // Current sequence position
        final int layer = 0;

        // Initialize input and output data
        FloatArray query = new FloatArray(dim);
        FloatArray keyCache = new FloatArray(seqLen * kvDim);
        FloatArray valueCache = new FloatArray(seqLen * kvDim);
        FloatArray output = new FloatArray(dim);
        FloatArray attentionWeights = new FloatArray(nHeads * (pos + 1));
        IntArray positionAndLayer = new IntArray(2); // Store position and layer indices

        // Initialize query vector with test values
        for (int i = 0; i < dim; i++) {
            query.set(i, 0.01f * i);
        }

        // Initialize key and value caches
        for (int i = 0; i < seqLen * kvDim; i++) {
            keyCache.set(i, 0.005f * i);
            valueCache.set(i, 0.005f * i);
        }

        // Clear output buffer
        for (int i = 0; i < dim; i++) {
            output.set(i, 0.0f);
        }

        // Set position and layer indices
        positionAndLayer.set(0, pos);
        positionAndLayer.set(1, layer);

        // Create reference implementation for comparison
        FloatArray expectedOutput = new FloatArray(dim);
        for (int i = 0; i < dim; i++) {
            expectedOutput.set(i, 0.0f);
        }

        // Execute reference CPU implementation
        processAttentionSequential(
                query, keyCache, valueCache, expectedOutput,
                nHeads, headSize, kvDim, kvMul, seqLen, pos, layer
        );

        // Create and execute TornadoVM task graph
        TaskGraph taskGraph = new TaskGraph("gqaTest")
                .transferToDevice(DataTransferMode.FIRST_EXECUTION,
                        query, keyCache, valueCache, output, attentionWeights, positionAndLayer)
                .task("parallel-attention", TornadoVMCompute::processHeadsParallel,
                        query, keyCache, valueCache, output,
                        nHeads, headSize, kvDim, kvMul, seqLen,
                        positionAndLayer, attentionWeights)
                .transferToHost(DataTransferMode.EVERY_EXECUTION, output);

        ImmutableTaskGraph immutableTaskGraph = taskGraph.snapshot();
        TornadoExecutionPlan executionPlan = new TornadoExecutionPlan(immutableTaskGraph);
        executionPlan.execute();

        // Verify results match the reference implementation
        float epsilon = 1e-4f;
        boolean resultsMatch = true;
        for (int i = 0; i < dim; i++) {
            if (Math.abs(output.get(i) - expectedOutput.get(i)) > epsilon) {
                resultsMatch = false;
                System.out.println("Mismatch at index " + i +
                        ": TornadoVM=" + output.get(i) +
                        ", Reference=" + expectedOutput.get(i));
            }
        }

        assertTrue("Parallel implementation results should match sequential implementation", resultsMatch);
    }
}