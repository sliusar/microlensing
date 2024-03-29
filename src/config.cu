#include <config.cuh>

Configuration::Configuration(const char* filename) {

    YAML::Node config = YAML::LoadFile(filename);

    verbose = config["verbose"].as<int>();
    
    sigma = config["sigma"].as<float>();
    sigma_c = config["sigma_c"].as<float>();
    gamma = config["gamma"].as<float>();
    R_field = config["R_field"].as<float>();
    M_max = config["M_max"].as<float>();
    V_max = config["V_max"].as<float>();
    
    R_rays = config["R_rays"].as<float>();
    dx_rays = config["dx_rays"].as<float>();
    
    dt = config["dt"].as<float>();
    t_max = config["t_max"].as<float>();
    
    image_width = config["image_width"].as<int>();
    image_height = config["image_height"].as<int>();
    
    image_y1_left = config["image_y1_left"].as<float>();
    image_y1_right = config["image_y1_right"].as<float>();
    
    image_y2_bottom = config["image_y2_bottom"].as<float>();
    image_y2_top = config["image_y2_top"].as<float>();

    save_images = config["save_images"].as<bool>();

    operation_mode = config["operation_mode"].as<int>();

    randomise_seed_number = config["randomise_seed_number"].as<int>();

    configuration_id = config["configuration_id"].as<string>();

    // Light curve
    lc_start_y1 = config["lc_start_y1"].as<float>();
    lc_start_y2 = config["lc_start_y2"].as<float>();
    lc_angle = config["lc_angle"].as<float>();
    lc_dist_max = config["lc_dist_max"].as<float>();
    lc_dist_step = config["lc_dist_step"].as<float>();

    lc_enabled = config["lc_enabled"].as<bool>();

    save_rays = config["save_rays"].as<bool>();

    std::vector<float> s = config["source_size"].as<std::vector<float>>();
    std::vector<float> e = config["eccentricity"].as<std::vector<float>>();
    
    std::copy(s.begin(), s.end(), source_size);
    std::copy(e.begin(), e.end(), eccentricity);

    p_ld = config["p_ld"].as<float>();
    p_pl = config["p_pl"].as<float>();

    // Recalculated values
    nMicrolenses = sigma * R_field * R_field / M_max;
    if (nMicrolenses == 0) nMicrolenses = 1;

    nRays = (int)ceil((M_PI * pow(R_rays, 2)) / pow(dx_rays, 2));
    nRaysLine = (int)ceil(2 * R_rays / dx_rays);
    nRaysSq = nRaysLine * nRaysLine;

    nCountSourceSizes = 0;
    nCountEccentricities = 0;
    nTimeSteps = 0;
    nLCsteps = 0;

    for (float s = source_size[0]; s < source_size[1]; s+=source_size[2]) nCountSourceSizes++;
    for (float e = eccentricity[0]; e < eccentricity[1]; e+=eccentricity[2]) nCountEccentricities++;
    for (float t = 0; t < t_max; t += dt) nTimeSteps++;
    for (float d = 0; d < lc_dist_max; d += lc_dist_step) nLCsteps++;

#if DEBUG == true 
    cout << endl << "DEBUGGING ACTIVATED: ADDITIONAL INFORMATION WILL BE INCLUDED IN THE LIGHT CURVE DATA" << endl << endl;
    nLCsteps += 1; // extra line for debugging information, see kernels.cu
#endif
    nLCcolumns = 2 + nCountSourceSizes * (2 /*ad*/ + 2 /*gs*/ + 2 /*ld*/ + 2 /*pl*/ + nCountEccentricities * (2 /*el*/ + 2 /*el_orth*/));

    image_pixel_y1_size = (image_y1_right - image_y1_left) / image_width;
    image_pixel_y2_size = (image_y2_top - image_y2_bottom) / image_height;
}

void Configuration::display() {
    cout << "CONFIGURATION ID: " << configuration_id << endl;
    cout << "--- GLS ---" << endl;
    cout << "  sigma: " << sigma << endl;
    cout << "  sigma_c: " << sigma_c << endl;
    cout << "  gamma: " << gamma << endl;
    cout << "  R_field: " << R_field << endl;
    cout << "  M_max: " << M_max << endl;
    cout << "  nMicrolenses: " << nMicrolenses << endl;
    cout << "  V_max: " << V_max << endl;
    cout << "  t_max: " << t_max << ", step: " << dt << endl;
    cout << endl;
  
    cout << "--- Ray tracing ---" << endl;
    cout << "  R_rays: " << R_rays << endl;
    cout << "  dx_rays: " << dx_rays << endl;
    cout << "  nRays: " << nRays << endl;
    cout << endl;
    
    if (lc_enabled) {
        cout << "--- LC trajectory ---" << endl;
        cout << "  Start y1: " << lc_start_y1 << endl;
        cout << "  Start y2: " << lc_start_y2 << endl;
        cout << "  Angle: " << lc_angle << endl;
        cout << "  T: [0, " << lc_dist_max << "], step: " << lc_dist_step << endl;
        cout << endl;
    }
    
    cout << "--- Image data ---" << endl;
    cout << "image_pixel_y1_size: " << image_pixel_y1_size << endl; 
    cout << "image_pixel_y2_size: " << image_pixel_y2_size << endl; 
    cout << "image_y1_left: " << image_y1_left << endl; 
    cout << "image_y1_right: " << image_y1_right << endl; 
    cout << "image_y2_bottom: " << image_y2_bottom << endl; 
    cout << "image_y2_top: " << image_y2_top << endl; 
    cout << endl;

    cout << "--- Sources ---" << endl;
    cout << "source_size: [" << source_size[0] << ", " << source_size[1] << "), step: " << source_size[2] << " (" << nCountSourceSizes << " steps in total)" << endl;
    cout << "eccentricity: [" << eccentricity[0] << ", " << eccentricity[1] << "), step: " << eccentricity[2] << " (" << nCountEccentricities << " steps in total)" << endl;
    cout << "p_ld: " << p_ld << endl;
    cout << "p_pl: " << p_pl << endl;

    cout << endl;
}

void Configuration::prepare_sources() {
    /*
    float src_size = source_size[0];
    float ecc = eccentricity[0];

    R_gs = src_size;
    R2_gs = R_gs * R_gs;

    p_ld = 2.0;
    R_1_2_ld = src_size * sqrt(log(2.0));
    R_ld = R_1_2_ld / sqrt(1.0 - pow(0.5, 2)/(p_ld + 2));
    R2_ld = R_ld * R_ld;

    p_pl = 3.0/2.0;
    R_1_2_pl = src_size * sqrt(log(2.0));
    R_pl = R_1_2_pl / sqrt((pow(2.0, 1.0/(p_pl - 1)) - 1.0)/log(2.0));
    R2_pl = R_pl * R_pl;

    R_1_2_ad = src_size * sqrt(log(2.0));
    R_ad = R_1_2_ad/4.0;
    R2_ad = R_ad * R_ad;

    e2_el = ecc * ecc;
    a_el = src_size / (1 - e2_el);
    b_el = src_size * (1 - e2_el);

    a2_el = a_el * a_el;
    b2_el = b_el * b_el;
    */
}