#ifndef COMMON_CUH
#define COMMON_CUH

#include <iostream>
#include <math.h>
#include <string>

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif 

#endif /* !COMMON_CUH */

// TODO: Separate C++ and Cuda, e.g. like in https://www.crossfire.nu/tutorials/179/mixing-cuda-and-c
