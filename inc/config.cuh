#ifndef INCLUDE_CONFIG_CUH
#define INCLUDE_CONFIG_CUH

#include <common.cuh>

class Configuration {
    char* filename;
    public:
        bool debug;

        bool save_images;
        int randomise_seed_number;
        float sigma, sigma_c, gamma, R_field, M_max, V_max, R_rays, dx_rays;
        float dt, t_max;
        int image_height, image_width;
        float image_y1_left, image_y2_bottom, image_y1_right, image_y2_top;

        string configuration_id;
        
        int nMicrolenses, nRays, nRaysSq, nRaysLine, nLCsteps, nTimeSteps; // Calculated        
        float image_pixel_y1_size, image_pixel_y2_size; // Calculated
        
        float lc_start_y1, lc_start_y2, lc_angle, lc_t_max, lc_t_step;
        bool lc_enabled, output_rays;

        float source_size;
        float R_gs, R2_gs, p_ld, R_1_2_ld, R_ld, R2_ld, p_pl, R_1_2_pl, R_pl, R2_pl, R_1_2_ad, R_ad, R2_ad; // Source-related, calculated

        Configuration(const char *);
        ~Configuration() {};
        void display();
        void prepare_sources();
};

#endif /* !INCLUDE_CONFIG_CUH */