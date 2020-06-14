#ifndef CONFIG_CUH
#define CONFIG_CUH

#include <common.cuh>
#include "yaml-cpp/yaml.h"

bool debug = true;

class Configuration {
    public:
        Configuration(char *);
        ~Configuration();
        void reconfigure();
};

#endif /* !CONFIG_CUH */