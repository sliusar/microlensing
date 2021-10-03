#ifndef INCLUDE_COMMON_CUH
#define INCLUDE_COMMON_CUH

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include "yaml-cpp/yaml.h"
#include <sys/types.h>
#include <sys/stat.h>

#ifdef __CUDACC__
    #define CUDA_CALLABLE_MEMBER __host__ __device__
#else
    #define CUDA_CALLABLE_MEMBER
#endif 

#define DEBUG false

struct __align__(16) Microlens { float x1, x2, v1, v2, m; } ;
struct __align__(16) Ray { float x1, x2; };
struct __align__(16) LC { float y1, y2, t, ampl; };

using namespace std;

#endif /* !INCLUDE_COMMON_CUH */

