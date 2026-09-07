#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <CL/cl.h>
#include <CL/cl_half.h>

// Define the matrix size
#define MATRIX_SIZE 4  

// Function to load kernel source code from a file
char* loadKernelSource(const char *filename, size_t *length) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Failed to load kernel file: %s\n", filename);
        exit(1);
    }
    fseek(file, 0, SEEK_END);
    *length = ftell(file);
    rewind(file);
    char *source = (char*)malloc(*length + 1);
    fread(source, 1, *length, file);
    source[*length] = '\0';
    fclose(file);
    return source;
}

int main() {
    cl_int err;

    // 1. Get platform and device
    cl_uint numPlatforms = 0;
    clGetPlatformIDs(0, NULL, &numPlatforms);  // Get the number of platforms
    cl_platform_id platforms[numPlatforms];
    clGetPlatformIDs(numPlatforms, platforms, NULL);  // Get platform list

    cl_platform_id selectedPlatform = NULL;
    cl_device_id selectedDevice = NULL;

    for (cl_uint i = 0; i < numPlatforms; i++) {
        char platformName[128];
        clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, sizeof(platformName), platformName, NULL);

        if (strstr(platformName, "Intel(R) OpenCL Graphics")) {
            selectedPlatform = platforms[i];

            // Try getting a GPU device from this platform
            err = clGetDeviceIDs(selectedPlatform, CL_DEVICE_TYPE_GPU, 1, &selectedDevice, NULL);
            if (err == CL_SUCCESS) {
                printf("Using platform: %s\n", platformName);
                break;
            }
        }
    }

    if (selectedPlatform == NULL || selectedDevice == NULL) {
        printf("Failed to find Intel GPU!\n");
        return 1;
    }

    // 2. Create OpenCL context and queue
    cl_context context = clCreateContext(NULL, 1, &selectedDevice, NULL, NULL, &err);
    if (err != CL_SUCCESS) {
        printf("Failed to create OpenCL context.\n");
        return 1;
    }

    cl_command_queue queue = clCreateCommandQueueWithProperties(context, selectedDevice, 0, &err);
    if (err != CL_SUCCESS) {
        printf("Failed to create OpenCL command queue.\n");
        return 1;
    }

    // 3. Create buffers
    size_t bufferSize = MATRIX_SIZE * sizeof(cl_half);
    cl_mem bufA = clCreateBuffer(context, CL_MEM_READ_ONLY, bufferSize, NULL, &err);
    cl_mem bufC = clCreateBuffer(context, CL_MEM_WRITE_ONLY, bufferSize, NULL, &err);
    if (!bufA || !bufC) {
        printf("Failed to allocate OpenCL buffers.\n");
        return 1;
    }

    // 4. Prepare input data (initialize with half values)
    cl_half inputA[MATRIX_SIZE];
    for (int i = 0; i < MATRIX_SIZE; i++) {
        inputA[i] = cl_half_from_float((float)(i + 1), CL_HALF_RTE);  // Convert float to half
    }

    // 5. Copy data to GPU
    err = clEnqueueWriteBuffer(queue, bufA, CL_TRUE, 0, bufferSize, inputA, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
        printf("Failed to write buffer to GPU.\n");
        return 1;
    }

    // 6. Load the kernel from file
    size_t kernelLength;
    char *kernelSource = loadKernelSource("floatconvert.cl", &kernelLength);

    // 7. Compile the kernel
    cl_program program = clCreateProgramWithSource(context, 1, (const char**)&kernelSource, &kernelLength, &err);
    err = clBuildProgram(program, 1, &selectedDevice, NULL, NULL, NULL);
    if (err != CL_SUCCESS) {
        char log[4096];
        clGetProgramBuildInfo(program, selectedDevice, CL_PROGRAM_BUILD_LOG, sizeof(log), log, NULL);
        printf("Kernel Build Error:\n%s\n", log);
        return 1;
    }

    // 8. Create kernel and set arguments
    cl_kernel kernel = clCreateKernel(program, "matrixMultiplicationHalfFloatToFloat", &err);
    if (err != CL_SUCCESS) {
        printf("Failed to create kernel.\n");
        return 1;
    }

    err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &bufA);
    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &bufC);
    if (err != CL_SUCCESS) {
        printf("Failed to set kernel arguments.\n");
        return 1;
    }

    // 9. Run the kernel
    size_t globalSize = MATRIX_SIZE;
    err = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &globalSize, NULL, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
        printf("Kernel execution failed.\n");
        return 1;
    }

    clFinish(queue);

    // 10. Read the result
    cl_half outputC[MATRIX_SIZE];
    err = clEnqueueReadBuffer(queue, bufC, CL_TRUE, 0, bufferSize, outputC, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
        printf("Failed to read buffer from GPU.\n");
        return 1;
    }

    // 11. Print results
    printf("Converted Half Precision values:\n");
    for (int i = 0; i < MATRIX_SIZE; i++) {
        printf("C[%d] = %f\n", i, cl_half_to_float(outputC[i]));
    }

    // Cleanup
    free(kernelSource);
    clReleaseMemObject(bufA);
    clReleaseMemObject(bufC);
    clReleaseKernel(kernel);
    clReleaseProgram(program);
    clReleaseCommandQueue(queue);
    clReleaseContext(context);

    return 0;
}
