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

#define CUDA_BLOCK_SIZE 512
#define CUDA_BLOCK_SIZE_2d 16

struct Microlens { float x1, x2, v1, v2, m; } ;
struct Ray { float x1, x2; };
struct Configuration { float sigma, sigma_c, gamma, R_field, M_avg, R_rays, dx_rays; int nMicrolenses, nRays; float dt, t_max; int image_height, image_width; };

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
    ul[i].x1 = x1;
    ul[i].x2 = x2;
    ul[i].v1 = 0.0;
    ul[i].v2 = 0.0;
    ul[i].m = 1.0;
  }
}

void populateRays(Ray *rays, int nRays, float R_rays, float dx_rays) {
  int counter = 0;
  for (float x1 = - R_rays; x1 <= R_rays; x1 += dx_rays) {
    for (float x2 = - R_rays; x2 <= R_rays; x2 += dx_rays) {
      rays[counter].x1 = x1;
      rays[counter].x2 = x2;
      counter += 1;
    }
  }
}

__device__ float distance2(float x, float y, float center_x, float center_y) {
  return pow(x - center_x, 2) + pow(y - center_y, 2);
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
      float m_x1 = ul.x1 + (ul.v1 * t);
      float m_x2 = ul.x2 + (ul.v2 * t);
      float rr = distance2(ray_x1, ray_x2, m_x1, m_x2);
      sum_x1 += ul.m * (ray_x1 - m_x1) / rr;
      sum_x2 += ul.m * (ray_x2 - m_x2) / rr;
    }
    rays[ri].x1 = (1 - c.gamma) * ray_x1 - c.sigma_c * ray_x1 - sum_x1;
    rays[ri].x2 = (1 + c.gamma) * ray_x2 - c.sigma_c * ray_x2 - sum_x2;
  }
}

__global__ void buildMap(Ray *rays, Configuration c) {
  int pix = blockDim.x * blockIdx.x + threadIdx.x;
  int piy = blockDim.y * blockIdx.y + threadIdx.y;
  if (pix < c.image_width && piy < c.image_height) {

  }
}

int main(const int argc, const char** argv) {
  srand (time(NULL));
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
  conf.nRays = pow((2 * conf.R_rays / conf.dx_rays) + 1, 2);

  // Iterations calculator
  conf.dt = 0.1;
  conf.t_max = 0.1;

  // Image definitions
  conf.image_height = 500;
  conf.image_width = 500;

  int ul_bytes = conf.nMicrolenses * sizeof(Microlens);
  int ray_bytes = conf.nRays * sizeof(Ray);

  Microlens *microlenses = (Microlens*)malloc(ul_bytes);
  Ray *rays = (Ray*)malloc(ray_bytes);

  randomiseMicrolenses(microlenses, conf.nMicrolenses, conf.R_field);
  populateRays(rays, conf.nRays, conf.R_rays, conf.dx_rays);

  Microlens *ul_buf;
  cudaMalloc(&ul_buf, ul_bytes);
  cudaMemcpy(ul_buf, microlenses, ul_bytes, cudaMemcpyHostToDevice);

  Ray *ray_buf;
  cudaMalloc(&ray_buf, ray_bytes);
  cudaMemcpy(ray_buf, rays, ray_bytes, cudaMemcpyHostToDevice);

  int nBlocks = (conf.nRays + CUDA_BLOCK_SIZE - 1) / CUDA_BLOCK_SIZE;

  ofstream f1, f2;
  f1.open("rays_x.dat");
  f2.open("rays_y.dat");
  for (float t = 0; t <= conf.t_max; t = t + conf.dt) {
    StartTimer();
    cout << "Time: " << t << endl;
    deflectRays<<<nBlocks, CUDA_BLOCK_SIZE>>>(ul_buf, ray_buf, conf, t); // compute ray deflections
    cudaMemcpy(rays, ray_buf, ray_bytes, cudaMemcpyDeviceToHost);
    cout << "Iteration completed in " << GetElapsedTime() << " ms" << endl;
    
    dim3 nBlocks_image = dim3((conf.image_width + CUDA_BLOCK_SIZE_2d - 1) / CUDA_BLOCK_SIZE_2d, (conf.image_height + CUDA_BLOCK_SIZE_2d - 1) / CUDA_BLOCK_SIZE_2d);
    cout << "Starting map calculation ..."<< endl;
    buildMap<<<nBlocks_image, dim3(CUDA_BLOCK_SIZE_2d, CUDA_BLOCK_SIZE_2d)>>>(ray_buf, conf); // build map from rays

    //cout << "Writing data ..."<< endl;
    //for (int i = 0; i <= conf.nRays; i++) f2 << rays[i].x1 << " " << rays[i].x2 << endl;
    //break;
  }
  f1.close();
  f2.close();

  free(microlenses);
  free(rays);
  
  cudaFree(ul_buf);
  cudaFree(ray_buf);
  return 0;
}
