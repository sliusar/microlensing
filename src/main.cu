#include <common.cuh>
#include <config.cuh>
#include <kernels.cuh>
#include "timer.h"

using namespace std;

#define cudaDeviceScheduleBlockingSync 0x04
#define cudaEventBlockingSync 0x01

#define CUDA_BLOCK_SIZE 1024
#define CUDA_BLOCK_SIZE_2d 32

#define LC_COLUMNS 4

int main(const int argc, const char** argv) {
  if (argc != 2) {
    cerr << "Usage:\n\t" << argv[0] << " configuration.yaml" << endl;
    return 1;
  }

  char filename[64], output_folder[64];

  Configuration conf(argv[1]);
  conf.display();
  if (conf.randomise_seed_number != 0) {
    long _seed = time(NULL);
    if (conf.randomise_seed_number > 0) _seed = conf.randomise_seed_number;
    cout << "Using " << _seed << " to seed the random generator" << endl;
    srand(_seed);
  }
  
  int ul_bytes = conf.nMicrolenses * sizeof(Microlens);
  int ray_bytes = conf.nRays * sizeof(Ray);
  int image_bytes = conf.image_height * conf.image_width * sizeof(int);
  int lc_bytes = LC_COLUMNS * conf.nLCsteps * sizeof(float);

  Microlens *microlenses = (Microlens*)malloc(ul_bytes);
  Ray *rays = (Ray*)malloc(ray_bytes);
  int *image = (int*)malloc(image_bytes);
  float *lc = (float*)malloc(lc_bytes);

  Microlens *ul_buf;
  Ray *ray_buf;
  int *image_buf;
  float *lc_buf;

  struct stat info;

  sprintf(output_folder, "./output/%s", conf.configuration_id.c_str());
  if( stat(output_folder , &info ) != 0 ) {
    if (mkdir(output_folder, 0755) != 0 && errno != EEXIST) {
      cerr << "Failed to create output folder " << output_folder << endl;
      return -1;
    }
  }

  cudaMalloc(&image_buf, image_bytes);
  
  cout << "Creating microlensing field ... " << flush;
  StartTimer();
  randomiseMicrolenses(microlenses, conf.nMicrolenses, conf.R_field);  
  cudaMalloc(&ul_buf, ul_bytes);
  cudaMemcpy(ul_buf, microlenses, ul_bytes, cudaMemcpyHostToDevice);
  cout << GetElapsedTime() << "s" << endl;

  cout << "Creating rays field ... " << flush;
  StartTimer();
  cudaMalloc(&ray_buf, ray_bytes);
  if (conf.output_rays) {
    populateRays(rays, conf.nRays, conf.R_rays, conf.dx_rays);
    cudaMemcpy(ray_buf, rays, ray_bytes, cudaMemcpyHostToDevice);
  }
  cout << GetElapsedTime() << "s" << endl;

  if (conf.lc_enabled) {
    cout << "Creating LC trajectory ... " << flush;
    StartTimer();
    createTrajectory(lc, conf);
    cudaMalloc(&lc_buf, lc_bytes);
    cudaMemcpy(lc_buf, lc, lc_bytes, cudaMemcpyHostToDevice);
    cout << GetElapsedTime() << "s" << endl;
  }

  int nBlocksRays = (conf.nRays + CUDA_BLOCK_SIZE - 1) / CUDA_BLOCK_SIZE;
  int nBlocksImageW = (conf.image_width + CUDA_BLOCK_SIZE_2d - 1) / CUDA_BLOCK_SIZE_2d;
  int nBlocksImageH = (conf.image_height + CUDA_BLOCK_SIZE_2d - 1) / CUDA_BLOCK_SIZE_2d;

  //From https://stackoverflow.com/questions/11888772/when-to-call-cudadevicesynchronize
  //  kernel1<<<X,Y>>>(...); // kernel start execution, CPU continues to next statement
  //  kernel2<<<X,Y>>>(...); // kernel is placed in queue and will start after kernel1 finishes, CPU continues to next statement
  //  cudaMemcpy(...); // CPU blocks until memory is copied, memory copy starts only after kernel2 finishes

  ofstream outf;
  int counter = 0;
  float _t = 0;
  float t_raytracing = 0;
  float t_output = 0;
  float t_lc = 0;
  for (float t = 0; t <= conf.t_max; t = t + conf.dt) {
    
    memset(image, 0, image_bytes);
    cudaMemcpy(image_buf, image, image_bytes, cudaMemcpyHostToDevice);

    cout << endl << ">>> Iteration #" << ++counter << ", t=" << t << endl;
    
    cout << "    Executing ray tracing ... " << flush;
    StartTimer();
    deflectRays<<<nBlocksRays, CUDA_BLOCK_SIZE>>>(ul_buf, ray_buf, conf, t, image_buf, lc_buf); // compute ray deflections
    if (conf.output_rays) cudaMemcpy(rays, ray_buf, ray_bytes, cudaMemcpyDeviceToHost);
    if (conf.lc_enabled) cudaMemcpy(lc, lc_buf, lc_bytes, cudaMemcpyDeviceToHost);
    if (conf.save_images) cudaMemcpy(image, image_buf, image_bytes, cudaMemcpyDeviceToHost);
    cudaError_t err = cudaDeviceSynchronize();
    if(err != cudaSuccess) {
      cerr << "Error running the deflectRays() kernel" << endl;
      return -1;
    }
    _t = GetElapsedTime();
    t_raytracing += _t;
    cout << _t << "s" << endl;

    
    if (conf.lc_enabled) {
      cout << "    Executing light curve calculation ... " << flush;
      StartTimer();
      calculateLCs<<<dim3(nBlocksImageW, nBlocksImageH), dim3(CUDA_BLOCK_SIZE_2d, CUDA_BLOCK_SIZE_2d)>>>(conf, image_buf, lc_buf); // compute lc
      cudaMemcpy(lc, lc_buf, lc_bytes, cudaMemcpyDeviceToHost);
      cudaError_t err = cudaDeviceSynchronize();
      if(err != cudaSuccess) {
        cerr << "Error running the calculateLCs() kernel" << endl;
        return -1;
      }
      _t = GetElapsedTime();
      t_lc += _t;
      cout << _t << "s" << endl;
    }

    //sprintf(filename, "%s/rays_y_%.2f.dat", output_folder, t);
    //cout << "  Writing data to " << filename << " ... " << flush;
    //outf.open(filename);
    //for (int i = 0; i <= conf.nRays; i++) {
    //  if (rays[i].d <= conf.image_diagonal_size && rays[i].x1 >= -30 && rays[i].x1 <= 30 && rays[i].x2 >= -30 && rays[i].x2 <= 30) {
    //    outf << rays[i].x1 << " " << rays[i].x2 << endl;
    //  }
    //}
    //outf.close();
    //cout << GetElapsedTime() << "s" << endl;

    sprintf(filename, "%s/image_%.2f.dat", output_folder, t);
    if (conf.save_images) {
      StartTimer();
      cout << "    Writing data to " << filename << " ... " << flush;
      outf.open(filename);
      outf << "# image (" << conf.image_width << ", " << conf.image_height << ")" << endl;
      outf << "# x in (" << conf.image_y1_left << ", " << conf.image_y1_right << ")" << endl;
      outf << "# y in (" << conf.image_y2_bottom << ", " << conf.image_y2_top << ")" << endl;

      for (int j = 0; j < conf.image_height; j++) {
        for (int i = 0; i < conf.image_width; i++) {
          outf << image[i * conf.image_width + j] << endl;
        }
      }
      outf.close();
      _t = GetElapsedTime();
      t_output += _t;
      cout << _t << "s" << endl;  
    } else {
      cout << "    Skipping writing data to " << filename << " ... " << endl;
    }

    sprintf(filename, "%s/lc_%.2f.dat", output_folder, t);
    if (conf.lc_enabled) {
      StartTimer();
      cout << "    Writing light curves data to " << filename << " ... " << flush;
      outf.open(filename);
      int c = conf.nLCsteps;
      for (int i = 0; i < c; i++) {
        outf << lc[i + 0 * c] << " " << lc[i + 1 * c] << " " << lc[i + 2 * c] << " " << lc[i + 3 * c] << endl;
      }
      outf.close();
      _t = GetElapsedTime();
      t_output += _t;
      cout << _t << "s" << endl;        
    }
  }

  cout << endl << ">>> Run summary:" << endl;
  cout << "    Raytracing time: " << t_raytracing << "s (mean: " << t_raytracing/counter << "s)" << endl;
  cout << "    Output time: " << t_output << "s (mean: " << t_output / counter << "s)" << endl;

  //if (conf.lc_enabled) printLC(lc, conf.nLCsteps);

  free(microlenses);
  free(rays);
  free(image);
  
  cudaFree(ul_buf);
  cudaFree(ray_buf);
  cudaFree(image_buf);

  return 0;
}
