#include <gravlens.hpp>


int main(int argc, char* argv[])
{
//  int devType=CL_DEVICE_TYPE_GPU;
//
//  if(argc > 1) {
//    devType = CL_DEVICE_TYPE_CPU;
//    Info("Using: CL_DEVICE_TYPE_CPU\n");
//  } else {
//    Info("Using: CL_DEVICE_TYPE_GPU\n");
//  }
    GravLens clGravLens(argc, argv);
    if(clGravLens.setupCL() != CL_SUCCESS) return -1;
    if(clGravLens.prepareData() != CL_SUCCESS) return -1;
    if(clGravLens.runCL() != CL_SUCCESS) return -1;

  return 0;
}

