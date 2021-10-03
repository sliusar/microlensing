#ifndef INCLUDE_KERNELS_CUH
#define INCLUDE_KERNELS_CUH

#include <common.cuh>

#include <common.cuh>
#include <config.cuh>

float distance(float, float, float, float);
float distance(float, float);
void randomiseMicrolenses(Microlens*, int, float);
void populateRays(Ray *, int, float, float);
void createTrajectory(float *, const Configuration);
void printLC(float *, int);

__device__ float dst2_inv(float, float);
__device__ float dst(float, float);
__device__ float dst(float, float, float, float);
__device__ float H(float, float);

__global__ void deflectRays(Microlens *, Ray *, const Configuration, const float, int *, float *);
__global__ void calculateLCs(const Configuration, int *, float *);

#endif /* !INCLUDE_KERNELS_CUH */