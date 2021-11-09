#include <common.cuh>
#include <config.cuh>
#include <kernels.cuh>
#include "timer.h"

using namespace std;

#define cudaDeviceScheduleBlockingSync 0x04
#define cudaEventBlockingSync 0x01

#define CUDA_BLOCK_SIZE 1024
#define CUDA_BLOCK_SIZE_2d 32

bool debug = false;

int write_image(char* filename, int* image, Configuration c) {
  FILE *fp = fopen(filename, "wb");
  if(fp == NULL) {
    cout << "Error opening the file " << filename << endl;
    return -1;
  }
  fwrite((const void*)&c.image_width, sizeof(c.image_width), 1, fp);
  fwrite((const void*)&c.image_height, sizeof(c.image_height), 1, fp);
  fwrite((const void*)&c.image_y1_left, sizeof(c.image_y1_left), 1, fp);
  fwrite((const void*)&c.image_y1_right, sizeof(c.image_y1_right), 1, fp);
  fwrite((const void*)&c.image_y2_bottom, sizeof(c.image_y2_bottom), 1, fp);
  fwrite((const void*)&c.image_y2_top, sizeof(c.image_y2_top), 1, fp);
  fwrite((const void*)&image[0], sizeof(image[0]), c.image_width * c.image_height, fp);
  fclose(fp);
  return 0;
}

double getCurrentTimestamp() {
  struct timeval time_now{};
  gettimeofday(&time_now, nullptr);
  return time_now.tv_sec + ((double)time_now.tv_usec / 1e6);
}

int estimateRaysCount(float R_rays, float dx_rays) {
  int counter = 0;
  for (float x1 = - R_rays; x1 <= R_rays; x1 += dx_rays) {
    for (float x2 = - R_rays; x2 <= R_rays; x2 += dx_rays) {
      if (distance(x1, x2) <= R_rays) counter++;
    }
  }
  return counter;
}

int main(const int argc, const char** argv) {
  double t0 = getCurrentTimestamp();
  if (argc != 2) {
    cerr << "Usage:\n\t" << argv[0] << " configuration.yaml" << endl;
    return 1;
  }

  char filename[64], output_folder[64];

  Configuration conf(argv[1]);
  conf.prepare_sources();
  conf.display();
  debug = conf.debug;
  if (conf.randomise_seed_number != 0) {
    long _seed = time(NULL);
    if (conf.randomise_seed_number > 0) _seed = conf.randomise_seed_number;
    cout << "Using " << _seed << " to seed the random generator" << endl;
    srand(_seed);
  }

  int _c = estimateRaysCount(conf.R_rays, conf.dx_rays);
  cout << "Print estimated rays count " << _c << " (previous " << conf.nRays << "). Adjusting value." << endl;
  conf.nRays = _c;

  int uls_bytes = conf.nMicrolenses * sizeof(Microlens);
  int rays_bytes = conf.nRays * sizeof(Ray);
  int image_bytes = conf.image_height * conf.image_width * sizeof(int);
  int lc_bytes = conf.nLCcolumns * conf.nLCsteps * sizeof(float);

  Microlens *microlenses = (Microlens*)malloc(uls_bytes);
  Ray *rays = (Ray*)malloc(rays_bytes);
  int *image = (int*)malloc(image_bytes);
  float *lc = (float*)malloc(lc_bytes);

  Microlens *ul_buf;
  Ray *rays_buf;
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
  
  cout << "Creating rays field ... " << flush;
  StartTimer();
  cudaMalloc(&rays_buf, rays_bytes);
  populateRays(rays, conf.nRays, conf.R_rays, conf.dx_rays);
  cudaMemcpy(rays_buf, rays, rays_bytes, cudaMemcpyHostToDevice);
  cout << GetElapsedTime() << "s" << endl;

  if (conf.lc_enabled) {
    cout << "Creating light curve placeholder ... " << flush;
    StartTimer();
    createLC(lc, conf);
    cudaMalloc(&lc_buf, lc_bytes);
    cudaMemcpy(lc_buf, lc, lc_bytes, cudaMemcpyHostToDevice);
    cout << GetElapsedTime() << "s" << endl;
  }

  int nBlocksRays = ceil((float)conf.nRays / (float)CUDA_BLOCK_SIZE);
  int nBlocksImageW = ceil((float)conf.image_width / (float)CUDA_BLOCK_SIZE_2d);
  int nBlocksImageH = ceil((float)conf.image_height / (float)CUDA_BLOCK_SIZE_2d);

  cout << endl << "GPU Execution settings:" << endl;
  cout << "    nBlocksRays: " << nBlocksRays << endl;
  cout << "    nBlocksImage: " << nBlocksImageW << " x " << nBlocksImageH << endl;

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
  bool microlenses_set = false;
  
  for (float t = 0; t <= conf.t_max; t = t + conf.dt) {
    cout << endl << ">>> Iteration #" << ++counter << ", t=" << t << " (elapsed: " << getCurrentTimestamp() - t0 << "s)" << endl;

    memset(image, 0, image_bytes);
    cudaMemcpy(image_buf, image, image_bytes, cudaMemcpyHostToDevice);

    resetLC(lc, conf);
    cudaMemcpy(lc_buf, lc, lc_bytes, cudaMemcpyHostToDevice);
    
    if (conf.operation_mode == 1 || (conf.operation_mode == 0 && microlenses_set == false)) {
      microlenses_set = true;
      cout << "    Creating microlensing field ... " << flush;
      StartTimer();
      randomiseMicrolenses(microlenses, conf);  
      cudaMalloc(&ul_buf, uls_bytes);
      cudaMemcpy(ul_buf, microlenses, uls_bytes, cudaMemcpyHostToDevice);
      cout << GetElapsedTime() << "s" << endl;
    }

    cout << "    [CUDA] Running ray tracing ... " << flush;
    StartTimer();
    deflectRays<<<nBlocksRays, CUDA_BLOCK_SIZE>>>(ul_buf, rays_buf, conf, t, image_buf, lc_buf); // compute ray deflections
    if (conf.output_rays) cudaMemcpy(rays, rays_buf, rays_bytes, cudaMemcpyDeviceToHost);
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
      cout << "    [CUDA] Calculating light curves ... " << flush;
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

    if (conf.output_rays) {
      sprintf(filename, "%s/rays_%.2f.dat", output_folder, t);
      cout << "    Writing data to " << filename << " ... " << flush;
      outf.open(filename);
      for (int i = 0; i <= conf.nRays; i++) {
        if (rays[i].x1 >= conf.image_y1_left && rays[i].x1 <= conf.image_y1_right && rays[i].x2 >= conf.image_y2_bottom && rays[i].x2 <= conf.image_y2_top) {
          outf << rays[i].x1 << " " << rays[i].x2 << endl;
        }
      }
      outf.close();
      cout << GetElapsedTime() << "s" << endl;
    }

    sprintf(filename, "%s/image_%.2f.dat", output_folder, t);
    if (conf.save_images) {
      StartTimer();
      cout << "    Writing data to " << filename << " ... " << flush;
      write_image(filename, image, conf);
      //outf.open(filename);
      //outf << "# image (" << conf.image_width << ", " << conf.image_height << ")" << endl;
      //outf << "# x in (" << conf.image_y1_left << ", " << conf.image_y1_right << ")" << endl;
      //outf << "# y in (" << conf.image_y2_bottom << ", " << conf.image_y2_top << ")" << endl;

      //for (int j = 0; j < conf.image_height; j++) {
      //  for (int i = 0; i < conf.image_width; i++) {
      //    outf << image[i * conf.image_width + j] << endl;
      //  }
      //}
      //outf.close();
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
 
      int counter = 0;
      outf << "# t y1 y2 ad gauss ld pl el el_orth" << endl;
      for (float t = 0.0; t < conf.lc_dist_max; t = t + conf.lc_dist_step) {
        if (counter < conf.nLCsteps) {
          outf << t << " ";
          outf << lc[counter + 0 * c] << " " << lc[counter + 1 * c] << " ";
          outf << lc[counter + 3 * c] / lc[counter + 2 * c] << " ";
          outf << lc[counter + 5 * c] / lc[counter + 4 * c] << " ";
          outf << lc[counter + 7 * c] / lc[counter + 6 * c] << " ";
          outf << lc[counter + 9 * c] / lc[counter + 8 * c] << " ";
          outf << lc[counter + 11 * c] / lc[counter + 10 * c] << " ";
          outf << lc[counter + 13 * c] / lc[counter + 12 * c] << " ";
          outf << endl;
        }
        counter++;
      }
      outf.close();
      _t = GetElapsedTime();
      t_output += _t;
      cout << _t << "s" << endl;        
    }
  }

  free(microlenses);
  free(rays);
  free(image);
  free(lc);
  
  cudaFree(ul_buf);
  cudaFree(lc_buf);
  cudaFree(rays_buf);
  cudaFree(image_buf);

  cout << endl << ">>> Run summary:" << endl;
  cout << "    Raytracing time: " << t_raytracing << "s (mean: " << t_raytracing/counter << "s)" << endl;
  cout << "    Output time: " << t_output << "s (mean: " << t_output / counter << "s)" << endl;
  cout << "    Total time: " << getCurrentTimestamp() - t0 << "s" << endl;

  return 0;
}
