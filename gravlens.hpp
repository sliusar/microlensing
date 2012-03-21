#include <iostream>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <CL/cl.h>

using namespace std;

#define __NO_STD_VECTOR // Use cl::vector and cl::string and
#define __NO_STD_STRING // not STL versions, more on this later

#define DATA_SIZE (1024*2098)

#define Info(...) fprintf(stdout, __VA_ARGS__)
#define Warn(...) fprintf(stderr, __VA_ARGS__)


const char* KernelFilename = "gravlens_kernel.cl";


static FILE * OpenTextFile(const char* filename, unsigned long* size /* returned file size in bytes */)
{
    struct stat statbuf;

    FILE *fh;
    fh = fopen(filename, "r");
    if (fh == 0)
        return NULL;

    stat(filename, &statbuf);
    if(size)
        (*size) = (unsigned long)statbuf.st_size;

    return fh;
}

static void CloseTextFile(FILE* fh)
{
    fclose(fh);
}


static unsigned long ReadFromTextFile(FILE* fh, char* buffer, size_t buffer_size)
{
    if(!fh)
        return 0;

    unsigned long count = (unsigned long)fread(buffer, buffer_size, 1, fh);
    buffer[buffer_size] = '\0';
    return count;
}

static char * LoadTextFromFile(const char *filename, unsigned long *size /* returned file size in bytes */)
{
    FILE* fh = OpenTextFile(filename, size);
    unsigned long bytes = (*size);

    char *text = (char*)malloc(bytes + 1);
    if(!text)
        return 0;

    ReadFromTextFile(fh, text, bytes);
    CloseTextFile(fh);

    return text;
}




class GravLens
{
    size_t                  global;
    size_t                   local;
    cl_int                     err;
    int                    devType;

    cl_platform_id      cpPlatform;
    cl_device_id         device_id;
    cl_context             context;
    cl_command_queue      commands;
    cl_program             program;
    cl_kernel               kernel;

    unsigned int           correct;
    float*                    data;
    float*                 results;

    cl_mem                   input;
    cl_mem                  output;

    unsigned int             count;

public:
    GravLens(int argc, char** argv) {
        devType = CL_DEVICE_TYPE_GPU;
        results = NULL;
        data = NULL;
    }
    ~GravLens() {
        delete [] data;
        delete [] results;

        clReleaseMemObject(input);
        clReleaseMemObject(output);
        clReleaseProgram(program);
        clReleaseKernel(kernel);
        clReleaseCommandQueue(commands);
        clReleaseContext(context);
    }

    int setupCL(void);
    int prepareData(void);
    int runCL(void);
};

int GravLens::setupCL() {
    // Connect to a compute device
    if (clGetPlatformIDs(1, &cpPlatform, NULL) != CL_SUCCESS) {
        Warn("Error: Failed to find a platform!\n");
        return EXIT_FAILURE;
    }

    // Get a device of the appropriate type
    if (clGetDeviceIDs(cpPlatform, devType, 1, &device_id, NULL) != CL_SUCCESS) {
        Warn("Error: Failed to create a device group!\n");
        return EXIT_FAILURE;
    }

    // Create a compute context
    context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
    if (!context) {
        Warn("Error: Failed to create a compute context!\n");
        return EXIT_FAILURE;
    }

    // Create a command commands
    commands = clCreateCommandQueue(context, device_id, 0, &err);
    if (!commands) {
        Warn("Error: Failed to create a command commands!\n");
        return EXIT_FAILURE;
    }

    // Load the kernel source code from disk
    //
    unsigned long kernel_length = 0;
    char* kernel_source = LoadTextFromFile(KernelFilename, &kernel_length);
    if (!kernel_source)
    {
        Warn("ERROR: Failed to load kernel from file!\n");
        return EXIT_FAILURE;
    }

    // Create the compute program from the source buffer
    program = clCreateProgramWithSource(context, 1, (const char **) &kernel_source, NULL, &err);
    if (!program || err != CL_SUCCESS) {
        Warn("Error: Failed to create compute program!\n");
        return EXIT_FAILURE;
    }

    // Build the program executable
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS) {
        size_t len;
        char buffer[2048];

        Warn("Error: Failed to build program executable!\n");
        clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        cerr << buffer << endl;
        exit(1);
    }

    // Create the compute kernel in the program
    kernel = clCreateKernel(program, "gravlens_kernel", &err);
    if (!kernel || err != CL_SUCCESS) {
        Warn("Error: Failed to create compute kernel!\n");
        exit(1);
    }

    return CL_SUCCESS;
}

int GravLens::prepareData() {
    count = DATA_SIZE;
    data = new float[DATA_SIZE];
    results = new float[DATA_SIZE];

    for(unsigned int i = 0; i < count; i++)
        data[i] = rand() / (float)RAND_MAX;

    // Create the device memory vectors
    //
    input = clCreateBuffer(context,  CL_MEM_READ_ONLY, sizeof(float) * count, NULL, NULL); 
    output = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * count, NULL, NULL);

    if (!input || !output) {
        Warn("Error: Failed to allocate device memory!\n");
        return EXIT_FAILURE;
    }

    // Transfer the input vector into device memory
    err = clEnqueueWriteBuffer(commands, input, CL_TRUE, 0, sizeof(float) * count, data, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
        Warn("Error: Failed to write to source array!\n");
        return EXIT_FAILURE;
    }
    return CL_SUCCESS;
}

int GravLens::runCL() {
    // Set the arguments to the compute kernel
    err = 0;
    err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &input);
    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &output);
    err |= clSetKernelArg(kernel, 2, sizeof(unsigned int), &count);

    if (err != CL_SUCCESS) {
        Warn("Error: Failed to set kernel arguments! %d", err);
        return EXIT_FAILURE;
    }

    // Get the maximum work group size for executing the kernel on the device
    err = clGetKernelWorkGroupInfo(kernel, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(local), &local, NULL);
    if (err != CL_SUCCESS) {
        Warn("Error: Failed to retrieve kernel work group info! %d", err);
        return EXIT_FAILURE;
    }

    Info("Work group size: %d\n", (int)local);
    // Execute the kernel over the vector using the 
    // maximum number of work group items for this device
    global = count;
    err = clEnqueueNDRangeKernel(commands, kernel, 1, NULL, &global, &local, 0, NULL, NULL);

    if (err) {
        Warn("Error: Failed to execute kernel!\n");
        return EXIT_FAILURE;
    }

    // Wait for all commands to complete
    clFinish(commands);

    // Read back the results from the device to verify the output
    //
    err = clEnqueueReadBuffer( commands, output, CL_TRUE, 0, sizeof(float) * count, results, 0, NULL, NULL );
    if (err != CL_SUCCESS) {
        Warn("Error: Failed to read output array! %d", err);
        exit(1);
    }

    // Validate our results
    //
    correct = 0;
//    for(unsigned int i = 0; i < count; i++) {
//        if(results[i] == data[i] * data[i])
//        correct++;
//    }

    // Print a brief summary detailing the results
//    cout << "Computed " << correct << "/" << count << " correct values" << endl;
//    cout << "Computed " << 100.f * (float)correct/(float)count << "% correct values" << endl;

    return CL_SUCCESS;
}
