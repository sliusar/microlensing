#ifndef INCLUDE_COMMON_CUH
#define INCLUDE_COMMON_CUH

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include "yaml-cpp/yaml.h"

#ifdef __CUDACC__
    #define CUDA_CALLABLE_MEMBER __host__ __device__
#else
    #define CUDA_CALLABLE_MEMBER
#endif 

#endif /* !INCLUDE_COMMON_CUH */
