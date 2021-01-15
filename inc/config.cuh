#ifndef CONFIG_CUH
#define CONFIG_CUH

#include <common.cuh>
#include "yaml-cpp/yaml.h"

class CudaConfiguration {
    int debug;
    public:
        CudaConfiguration(char *);
        ~CudaConfiguration();
        void reconfigure();
};

#endif /* !CONFIG_CUH */