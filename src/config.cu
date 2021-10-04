#include <config.cuh>

Configuration::Configuration(const char* filename) {

    YAML::Node config = YAML::LoadFile(filename);
    
    sigma = config["sigma"].as<float>();
    sigma_c = config["sigma_c"].as<float>();
    gamma = config["gamma"].as<float>();
    R_field = config["R_field"].as<float>();
    M_avg = config["M_avg"].as<float>();
    
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
    lc_end_y1 = config["lc_end_y1"].as<float>();
    lc_end_y2 = config["lc_end_y2"].as<float>();
    lc_step = config["lc_step"].as<float>();

    lc_enabled = config["lc_enabled"].as<bool>();

    output_rays = config["output_rays"].as<bool>();

    // Recalculated values
    nMicrolenses = sigma * R_field * R_field / M_avg;
    if (nMicrolenses == 0) nMicrolenses = 1;

    nRays = (int)ceil((M_PI * pow(R_rays, 2)) / pow(dx_rays, 2));
    nRaysLine = (int)ceil(2 * R_rays / dx_rays);
    nRaysSq = nRaysLine * nRaysLine;

    nLCsteps = (int)ceil(sqrt(pow(lc_end_y1 - lc_start_y1, 2) + pow(lc_end_y2 - lc_start_y2, 2))/lc_step);
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
    cout << "  M_avg: " << M_avg << endl;
    cout << "  nMicrolenses: " << nMicrolenses << endl;
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
        cout << "  End y1: " << lc_end_y1 << endl;
        cout << "  End y2: " << lc_end_y2 << endl;
        cout << "  Step: " << lc_step << endl;
        cout << "  LC steps: " << nLCsteps << endl;    
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
}