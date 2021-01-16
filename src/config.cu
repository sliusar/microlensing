#include <config.cuh>
using namespace std;

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
    
    image_y_width = config["image_y_width"].as<float>();
    image_y_height = config["image_y_height"].as<float>();
    
    image_center_y1 = config["image_center_y1"].as<float>();
    image_center_y2 = config["image_center_y2"].as<float>();

    randomise_with_time = config["randomise_with_time"].as<bool>();

    // Recalculated values
    nMicrolenses = sigma * M_PI * R_field * R_field / M_PI * M_avg;

    nRays = (int)ceil((M_PI * pow(R_rays, 2)) / pow(dx_rays, 2));

    image_pixel_y1_size = image_y_width / image_width;
    image_pixel_y2_size = image_y_height / image_height;
  
    image_y1_left = image_center_y1 - image_y_width/2;
    image_y1_right = image_center_y1 + image_y_width/2;
    image_y2_bottom = image_center_y2 - image_y_height/2;
    image_y2_top = image_center_y2 + image_y_height/2;
}

void Configuration::display() {
    cout << "CONFIGURATION:" << endl;
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
}