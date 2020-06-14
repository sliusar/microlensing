#ifndef CONFIG_CUH
#define CONFIG_CUH

#include <common.cuh>
#include "yaml-cpp/yaml.h"

class Configuration {
    int debug;
    public:
        Configuration(char *);
        ~Configuration();
        void reconfigure();
};

#endif /* !CONFIG_CUH */