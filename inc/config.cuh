#ifndef INCLUDE_CONFIG_CUH
#define INCLUDE_CONFIG_CUH

#include <common.cuh>
using namespace std;

class Configuration {
    int debug;
    char* filename;
    public:
        bool randomise_with_time;
        float sigma, sigma_c, gamma, R_field, M_avg, R_rays, dx_rays;
        float dt, t_max;
        int image_height, image_width;
        float image_y2_height, image_y1_width, image_center_y1, image_center_y2;

        string configuration_id;
        
        int nMicrolenses, nRays, nRays_square, nRays_line; // Calculated        
        float image_pixel_y1_size, image_pixel_y2_size, image_y1_left, image_y2_bottom, image_y1_right, image_y2_top, image_diagonal_size; // Calculated
        
        Configuration(const char *);
        ~Configuration() {};
        void reconfigure();
        void setdebug(bool);
        void display();
};

#endif /* !INCLUDE_CONFIG_CUH */