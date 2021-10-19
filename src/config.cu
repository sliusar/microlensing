#include <config.cuh>

Configuration::Configuration(const char* filename) {

    YAML::Node config = YAML::LoadFile(filename);

    debug = config["debug"].as<bool>();
    
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

    randomise_seed_number = config["randomise_seed_number"].as<int>();

    configuration_id = config["configuration_id"].as<string>();

    // Light curve
    lc_start_y1 = config["lc_start_y1"].as<float>();
    lc_start_y2 = config["lc_start_y2"].as<float>();
    lc_angle = config["lc_angle"].as<float>();
    lc_t_max = config["lc_t_max"].as<float>();
    lc_t_step = config["lc_t_step"].as<float>();

    lc_enabled = config["lc_enabled"].as<bool>();

    output_rays = config["output_rays"].as<bool>();

    source_size = config["source_size"].as<float>();

    // Recalculated values
    nMicrolenses = sigma * R_field * R_field / M_max;
    if (nMicrolenses == 0) nMicrolenses = 1;

    nRays = (int)ceil((M_PI * pow(R_rays, 2)) / pow(dx_rays, 2));
    nRaysLine = (int)ceil(2 * R_rays / dx_rays);
    nRaysSq = nRaysLine * nRaysLine;

    nLCsteps = (int)round(lc_t_max / lc_t_step);
    nTimeSteps = (int)round((t_max + dt) / dt);

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
    cout << "  t_max: " << t_max << endl;
    cout << "  dt: " << dt << endl;
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
        cout << "  T: [0, " << lc_t_max << "]" << endl;
        cout << "  dT: " << lc_t_step << endl;  
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
    cout << "source_size: " << source_size << endl;
    cout << "R_gs: " << R_gs << endl;
    cout << "p_ld: " << p_ld << endl;
    cout << "R_1_2_ld: " << R_1_2_ld << endl;
    cout << "R_ld: " << R_ld << endl;
    cout << "p_pl: " << p_pl << endl;
    cout << "R_1_2_pl: " << R_1_2_pl << endl;
    cout << "R_pl: " << R_pl << endl;
    cout << "R_1_2_ad: " << R_1_2_ad << endl;
    cout << "R_ad: " << R_ad << endl;
    cout << endl;
}

void Configuration::prepare_sources() {

    R_gs = source_size;
    R2_gs = R_gs * R_gs;

    p_ld = 2.0;
    R_1_2_ld = source_size * sqrt(log(2.0));
    R_ld = R_1_2_ld / sqrt(1.0 - pow(0.5, 2)/(p_ld + 2));
    R2_ld = R_ld * R_ld;

    p_pl = 3.0/2.0;
    R_1_2_pl = source_size * sqrt(log(2.0));
    R_pl = R_1_2_pl / sqrt((pow(2.0, 1.0/(p_pl - 1)) - 1.0)/log(2.0));
    R2_pl = R_pl * R_pl;

    R_1_2_ad = source_size * sqrt(log(2.0));
    R_ad = R_1_2_ad/4.0;
    R2_ad = R_ad * R_ad;
}