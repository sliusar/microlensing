#include <common.cuh>
#include <config.cuh>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "timer.h"
#include <time.h>

#include <iostream>
#include <fstream>
using namespace std;

#define CUDA_BLOCK_SIZE 1024
#define CUDA_BLOCK_SIZE_2d 32

struct Microlens { float x1, x2, v1, v2, m; } ;
struct Ray { float x1, x2; };
struct Configuration { float sigma, sigma_c, gamma, R_field, M_avg, R_rays, dx_rays; int nMicrolenses, nRays, nRays1d; float dt, t_max; int image_height, image_width; float image_y_height, image_y_width, image_center_y1, image_center_y2, image_pixel_y1_size, image_pixel_y2_size, image_y1_left, image_y2_bottom, image_y1_right, image_y2_top; };

float distance(float x, float y, float center_x, float center_y) {
  return sqrt(pow(x - center_x, 2) + pow(y - center_y, 2));
}

void randomiseMicrolenses(Microlens *ul, int n, float R) {
  for (int i = 0; i < n; i++) {
    float x1 = 2 * R * (rand() / (float)RAND_MAX) - R;
    float x2 = 2 * R * (rand() / (float)RAND_MAX) - R;
    while (distance(x1, x2, 0, 0) > R) {
      x1 = 2 * R * (rand() / (float)RAND_MAX) - R;
      x2 = 2 * R * (rand() / (float)RAND_MAX) - R;
    }
    ul[i] = {.x1 = x1, .x2 = x2, .v1 = 0.0, .v2 = 0.0, .m = 1.0 };
  }
}

void populateRays(Ray *rays, int nRays, float R_rays, float dx_rays) {
  int counter = 0;
  for (float x1 = - R_rays; x1 <= R_rays; x1 += dx_rays) {
    for (float x2 = - R_rays; x2 <= R_rays; x2 += dx_rays) {
      rays[counter++] = {.x1 = x1, .x2 = x2};
    }
  }
}

__device__ float dst2_inv(float x, float y) {
  return pow(rhypotf(x, y), 2);
}

__global__ void deflectRays(Microlens *uls, Ray *rays, const Configuration c, const float t) {
  int ri = blockDim.x * blockIdx.x + threadIdx.x;
  if (ri < c.nRays) {
    float ray_x1 = rays[ri].x1;
    float ray_x2 = rays[ri].x2;
    float sum_x1 = 0.0;
    float sum_x2 = 0.0;
    for (int i = 0; i < c.nMicrolenses; i++) {
      Microlens ul = uls[i];
      float m_x1 = ray_x1 - ul.x1 - (ul.v1 * t);
      float m_x2 = ray_x2 - ul.x2 - (ul.v2 * t);
      float ri = ul.m * dst2_inv(m_x1, m_x2);
      sum_x1 += m_x1 * ri;
      sum_x2 += m_x2 * ri;
    }
    rays[ri].x1 = (1 - c.gamma) * ray_x1 - c.sigma_c * ray_x1 - sum_x1;
    rays[ri].x2 = (1 + c.gamma) * ray_x2 - c.sigma_c * ray_x2 - sum_x2;
  }
}

__global__ void buildMap(Ray *rays, const Configuration c, float *image) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  int j = blockDim.y * blockIdx.y + threadIdx.y;
  if (i < c.image_width && j < c.image_height) {
    for (int r = 0; i < c.nRays; i++) {

    }
  }
}

int main(const int argc, const char** argv) {
  //srand (time(NULL));
  Configuration conf;
  // Microlensing field configuration
  conf.sigma = 0.5;
  conf.sigma_c = 0;
  conf.gamma = 0.1;
  conf.R_field = 100;
  conf.M_avg = 1.0;
  conf.nMicrolenses = conf.sigma * M_PI * conf.R_field * conf.R_field / M_PI * conf.M_avg;
  printf("sigma: %f\nsigma_c: %f\ngamma: %f\nR: %f\nM_avg: %f\nN: %d\n", conf.sigma, conf.sigma_c, conf.gamma, conf.R_field, conf.M_avg, conf.nMicrolenses);

  // Rays configuration
  conf.R_rays = 100;
  conf.dx_rays = 0.1;
  conf.nRays1d = 2 * conf.R_rays / conf.dx_rays + 1;
  conf.nRays = pow(conf.nRays1d, 2);

  // Iterations calculator
  conf.dt = 0.1;
  conf.t_max = 0.1;

  // Image definitions
  conf.image_width = 500;
  conf.image_height = 500;
  
  conf.image_y_width = 30.0;
  conf.image_y_height = 30.0;
  
  conf.image_center_y1 = 0.0;
  conf.image_center_y2 = 0.0;

  conf.image_pixel_y1_size = conf.image_y_width / conf.image_width;
  conf.image_pixel_y2_size = conf.image_y_height / conf.image_height;

  conf.image_y1_left = conf.image_center_y1 - conf.image_y_width/2;
  conf.image_y1_right = conf.image_center_y1 + conf.image_y_width/2;
  conf.image_y2_bottom = conf.image_center_y2 - conf.image_y_height/2;
  conf.image_y2_top = conf.image_center_y2 + conf.image_y_height/2;

  int ul_bytes = conf.nMicrolenses * sizeof(Microlens);
  int ray_bytes = conf.nRays * sizeof(Ray);
  int image_bytes = conf.image_height * conf.image_width * sizeof(float);

  Microlens *microlenses = (Microlens*)malloc(ul_bytes);
  Ray *rays = (Ray*)malloc(ray_bytes);
  //float *image = (float*)malloc(image_bytes);

  randomiseMicrolenses(microlenses, conf.nMicrolenses, conf.R_field);
  populateRays(rays, conf.nRays, conf.R_rays, conf.dx_rays);
  //memset(image, 0, sizeof(image));

  Microlens *ul_buf;
  cudaMalloc(&ul_buf, ul_bytes);
  cudaMemcpy(ul_buf, microlenses, ul_bytes, cudaMemcpyHostToDevice);

  Ray *ray_buf;
  cudaMalloc(&ray_buf, ray_bytes);
  cudaMemcpy(ray_buf, rays, ray_bytes, cudaMemcpyHostToDevice);

  //float *image_buf;
  //cudaMalloc(&image_buf, image_bytes);
  //cudaMemcpy(image_buf, image, image_bytes, cudaMemcpyHostToDevice);

  int nBlocks = (conf.nRays + CUDA_BLOCK_SIZE - 1) / CUDA_BLOCK_SIZE;

  ofstream outf;
  for (float t = 0; t <= conf.t_max; t = t + conf.dt) {
    StartTimer();
    deflectRays<<<nBlocks, CUDA_BLOCK_SIZE>>>(ul_buf, ray_buf, conf, t); // compute ray deflections
    cudaMemcpy(rays, ray_buf, ray_bytes, cudaMemcpyDeviceToHost);
    cout << "Iteration t=" << t << " completed in " << GetElapsedTime() << " ms" << endl;
    
    cudaDeviceSynchronize();
    dim3 nBlocks_image = dim3((conf.image_width + CUDA_BLOCK_SIZE_2d - 1) / CUDA_BLOCK_SIZE_2d, (conf.image_height + CUDA_BLOCK_SIZE_2d - 1) / CUDA_BLOCK_SIZE_2d);
    cout << "Starting map calculation ..."<< endl;
    //buildMap<<<nBlocks_image, dim3(CUDA_BLOCK_SIZE_2d, CUDA_BLOCK_SIZE_2d)>>>(ray_buf, conf, image_buf); // build map from rays

    char filename[32];
    sprintf(filename, "rays_y_%.2f.dat", t);
    cout << "Writing data to " << filename << " ..."<< endl;
    outf.open(filename);
    for (int i = 0; i <= conf.nRays; i++) {
      if (rays[i].x1 >= -20 && rays[i].x1 <= 20) {
        if (rays[i].x2 >= -20 && rays[i].x2 <= 20) {
          outf << rays[i].x1 << " " << rays[i].x2 << endl;
        }
      }
    }
    outf.close();
    break;
  }

  free(microlenses);
  free(rays);
  
  cudaFree(ul_buf);
  cudaFree(ray_buf);
  return 0;
}
