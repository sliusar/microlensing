#include <common.cuh>
#include <config.cuh>
#include "timer.h"

using namespace std;

#define CUDA_BLOCK_SIZE 1024
#define CUDA_BLOCK_SIZE_2d 32

struct Microlens { float x1, x2, v1, v2, m; } ;
struct Ray { float x1, x2; };

float distance(float x, float y, float center_x, float center_y) {
  return sqrt(pow(x - center_x, 2) + pow(y - center_y, 2));
}

float distance(float x, float y) {
  return distance(x, y, 0, 0);
}

void randomiseMicrolenses(Microlens *ul, int n, float R) {
  for (int i = 0; i < n; i++) {
    float x1 = 2 * R * (rand() / (float)RAND_MAX) - R;
    float x2 = 2 * R * (rand() / (float)RAND_MAX) - R;
    while (distance(x1, x2) > R) {
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
      if (distance(x1, x2) <= R_rays && counter < nRays) rays[counter++] = {.x1 = x1, .x2 = x2};
    }
  }
}

__device__ float dst2_inv(float x, float y) {
  return pow(rhypotf(x, y), 2);
}

float dst2_inv_cpu(float x, float y) {
  return 1 / (pow(x, 2) + pow(y, 2));
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

void deflectRaysCPU(Microlens *uls, Ray *rays, const Configuration c, const float t) {
  for (int ri = 0; ri < c.nRays; ri++) {
    float ray_x1 = rays[ri].x1;
    float ray_x2 = rays[ri].x2;
    float sum_x1 = 0.0;
    float sum_x2 = 0.0;
    for (int i = 0; i < c.nMicrolenses; i++) {
      Microlens ul = uls[i];
      float m_x1 = ray_x1 - ul.x1 - (ul.v1 * t);
      float m_x2 = ray_x2 - ul.x2 - (ul.v2 * t);
      float ri = ul.m * dst2_inv_cpu(m_x1, m_x2);
      sum_x1 += m_x1 * ri;
      sum_x2 += m_x2 * ri;
    }
    rays[ri].x1 = (1 - c.gamma) * ray_x1 - c.sigma_c * ray_x1 - sum_x1;
    rays[ri].x2 = (1 + c.gamma) * ray_x2 - c.sigma_c * ray_x2 - sum_x2;
  }
}

__global__ void buildMap(Ray *rays, const Configuration c, float *image) {
  int row = blockIdx.y * blockDim.y + threadIdx.y;
  int col = blockIdx.x * blockDim.x + threadIdx.x;
  if (col < c.image_width && row < c.image_height) {
    float y1_min = c.image_y1_left + col * c.image_pixel_y1_size;
    float y1_max = y1_min + c.image_pixel_y1_size;
    float y2_min = c.image_y2_bottom + row * c.image_pixel_y2_size;
    float y2_max = y2_min + c.image_pixel_y1_size;
    for (int i = 0; i < c.nRays; i++) {
      Ray ray = rays[i];
      if (y1_min <= ray.x1 && ray.x1 < y1_max) {
        if (y2_min <= ray.x2 && ray.x2 < y2_max) {
          image[col * c.image_width + row]++;
        }
      }
    }
  }
}

int main(const int argc, const char** argv) {
  if (argc != 2) {
    cerr << "Usage:\n\t" << argv[0] << " configuration.yaml" << endl;
    return 1;
  }

  Configuration conf(argv[1]);
  conf.display();
  if (conf.randomise_with_time) srand(time(NULL));
  
  int ul_bytes = conf.nMicrolenses * sizeof(Microlens);
  int ray_bytes = conf.nRays * sizeof(Ray);
  int image_bytes = conf.image_height * conf.image_width * sizeof(float);

  Microlens *microlenses = (Microlens*)malloc(ul_bytes);
  Ray *rays = (Ray*)malloc(ray_bytes);
  float *image = (float*)malloc(image_bytes);

  StartTimer();
  randomiseMicrolenses(microlenses, conf.nMicrolenses, conf.R_field);  
  cout << "Creating microlensing field in " << GetElapsedTime() << " s" << endl;
  populateRays(rays, conf.nRays, conf.R_rays, conf.dx_rays);
  cout << "Defining rays field in " << GetElapsedTime() << " s" << endl;
  memset(image, 0, sizeof(image));

  Microlens *ul_buf;
  cudaMalloc(&ul_buf, ul_bytes);
  cudaMemcpy(ul_buf, microlenses, ul_bytes, cudaMemcpyHostToDevice);

  Ray *ray_buf;
  cudaMalloc(&ray_buf, ray_bytes);
  cudaMemcpy(ray_buf, rays, ray_bytes, cudaMemcpyHostToDevice);

  float *image_buf;
  cudaMalloc(&image_buf, image_bytes);
  cudaMemcpy(image_buf, image, image_bytes, cudaMemcpyHostToDevice);

  int nBlocks = (conf.nRays + CUDA_BLOCK_SIZE - 1) / CUDA_BLOCK_SIZE;

  ofstream outf;
  for (float t = 0; t <= conf.t_max; t = t + conf.dt) {
    cout << "Iteration t: " << t << endl;
    cout << "  Executing ray tracing ... " << endl;
    StartTimer();
    deflectRays<<<nBlocks, CUDA_BLOCK_SIZE>>>(ul_buf, ray_buf, conf, t); // compute ray deflections
    //deflectRaysCPU(microlenses, rays, conf, t); // CPU version
    cudaMemcpy(rays, ray_buf, ray_bytes, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cout << "    done in " << GetElapsedTime() << " s" << endl;
    

    cout << "  Executing amplification map calculation ... " << endl;
    StartTimer();
    dim3 nBlocks_image = dim3((conf.image_width + CUDA_BLOCK_SIZE_2d - 1) / CUDA_BLOCK_SIZE_2d, (conf.image_height + CUDA_BLOCK_SIZE_2d - 1) / CUDA_BLOCK_SIZE_2d);
    buildMap<<<nBlocks_image, dim3(CUDA_BLOCK_SIZE_2d, CUDA_BLOCK_SIZE_2d)>>>(ray_buf, conf, image_buf); // build map from rays
    cudaMemcpy(image, image_buf, image_bytes, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cout << "    done in " << GetElapsedTime() << " s" << endl;

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

    outf.open("image.dat");
    for (int j = 0; j < conf.image_height; j++) {
      for (int i = 0; i < conf.image_width; i++) {
        outf << image[i * conf.image_width + j] << endl;
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
