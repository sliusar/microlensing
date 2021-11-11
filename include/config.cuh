#ifndef INCLUDE_CONFIG_CUH
#define INCLUDE_CONFIG_CUH

#include <common.cuh>

class Configuration {
    char* filename;
    public:
        int verbose;

        int operation_mode;
        bool save_images;
        int randomise_seed_number;
        float sigma, sigma_c, gamma, R_field, M_max, V_max, R_rays, dx_rays;
        float dt, t_max;
        int image_height, image_width;
        float image_y1_left, image_y2_bottom, image_y1_right, image_y2_top;

        string configuration_id;
        
        int nMicrolenses, nRays, nRaysSq, nRaysLine, nLCsteps, nLCcolumns, nTimeSteps; // Calculated        
        float image_pixel_y1_size, image_pixel_y2_size; // Calculated
        
        float lc_start_y1, lc_start_y2, lc_angle, lc_dist_max, lc_dist_step;
        bool lc_enabled, save_rays;

        float source_size[3];
        float eccentricity[3];
        int nCountSourceSizes;
        int nCountEccentricities;

        float p_ld, p_pl;

        //float R_gs, R2_gs, R_1_2_ld, R_ld, R2_ld, R_1_2_pl, R_pl, R2_pl, R_1_2_ad, R_ad, R2_ad; // Source-related, calculated
        //float e_el, e2_el, a_el, b_el, a2_el, b2_el;

        Configuration(const char *);
        ~Configuration() {};
        void display();
        void prepare_sources();
};

#endif /* !INCLUDE_CONFIG_CUH */