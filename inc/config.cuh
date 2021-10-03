#ifndef INCLUDE_CONFIG_CUH
#define INCLUDE_CONFIG_CUH

#include <common.cuh>

class Configuration {
    int debug;
    char* filename;
    public:
        bool save_images;
        int randomise_seed_number;
        float sigma, sigma_c, gamma, R_field, M_avg, R_rays, dx_rays;
        float dt, t_max;
        int image_height, image_width;
        float image_y1_left, image_y2_bottom, image_y1_right, image_y2_top;

        string configuration_id;
        
        int nMicrolenses, nRays, nRays_square, nRays_line, nLCsteps, nTimeSteps; // Calculated        
        float image_pixel_y1_size, image_pixel_y2_size; // Calculated
        
        float lc_start_y1, lc_start_y2, lc_end_y1, lc_end_y2, lc_step;
        bool lc_enabled, output_rays;

        Configuration(const char *);
        ~Configuration() {};
        void reconfigure();
        void setdebug(bool);
        void display();
};

#endif /* !INCLUDE_CONFIG_CUH */