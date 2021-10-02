#include <common.cuh>
#include <config.cuh>
#include "timer.h"

using namespace std;

#define cudaDeviceScheduleBlockingSync 0x04
#define cudaEventBlockingSync 0x01

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

  float speed_range_radius = 1;
  for (int i = 0; i < n; i++) {
    float v1 = speed_range_radius * (rand() / (float)RAND_MAX) - speed_range_radius;
    float v2 = speed_range_radius * (rand() / (float)RAND_MAX) - speed_range_radius;
    ul[i].v1 = v1;
    ul[i].v2 = v2;
  }
}

void populateRays(Ray *rays, int nRays, float R_rays, float dx_rays) {
  int counter = 0;
  for (float x1 = - R_rays; x1 <= R_rays; x1 += dx_rays) {
    for (float x2 = - R_rays; x2 <= R_rays; x2 += dx_rays) {
      if (distance(x1, x2) <= R_rays && counter < nRays) rays[counter++] = {.x1 = x1, .x2 = x2 };
    }
  }
}

__device__ float dst2_inv(float x, float y) {
  return rhypotf(x, y) * rhypotf(x, y);
}

__global__ void createRays(Ray *rays, const Configuration c) {
  int ri = blockDim.x * blockIdx.x + threadIdx.x;
  if (ri < c.nRays_square) {
    int j = ri / c.nRays_line;
    int i = ri - j * c.nRays_line;

    rays[ri].x1 = i * c.dx_rays - c.R_rays;
    rays[ri].x2 = j * c.dx_rays - c.R_rays;
  }
}

__global__ void deflectRays(Microlens *uls, Ray *rays, const Configuration c, const float t, int *image) {
  int ri = blockDim.x * blockIdx.x + threadIdx.x;
  if (ri < c.nRays) {
    float ray_x1 = rays[ri].x1;
    float ray_x2 = rays[ri].x2;
    float sum_x1 = 0.0;
    float sum_x2 = 0.0;
    for (int i = 0; i < c.nMicrolenses; i++) {
      float m_x1 = ray_x1 - uls[i].x1 - (uls[i].v1 * t);
      float m_x2 = ray_x2 - uls[i].x2 - (uls[i].v2 * t);
      float ri = uls[i].m * dst2_inv(m_x1, m_x2);
      sum_x1 += m_x1 * ri;
      sum_x2 += m_x2 * ri;
    }
    //rays[ri].x1 = (1 - c.gamma) * ray_x1 - c.sigma_c * ray_x1 - sum_x1;
    //rays[ri].x2 = (1 + c.gamma) * ray_x2 - c.sigma_c * ray_x2 - sum_x2;
    //rays[ri].d = hypotf(rays[ri].x1 - c.image_center_y1, rays[ri].x2 - c.image_center_y2);
    //int x = lrintf((rays[ri].x1 - c.image_y1_left) / c.image_pixel_y1_size);
    //int y = lrintf((rays[ri].x2 - c.image_y2_bottom) / c.image_pixel_y2_size);

    float r_x1 = (1 - c.gamma) * ray_x1 - c.sigma_c * ray_x1 - sum_x1;
    float r_x2 = (1 + c.gamma) * ray_x2 - c.sigma_c * ray_x2 - sum_x2;
    int x = lrintf((r_x1 - c.image_y1_left) / c.image_pixel_y1_size);
    int y = lrintf((r_x2 - c.image_y2_bottom) / c.image_pixel_y2_size);
    if (x >= 0 && x < c.image_width && y >= 0 && y < c.image_height) atomicAdd(&image[x * c.image_width + y], 1.0);
  }
}

int main(const int argc, const char** argv) {
  if (argc != 2) {
    cerr << "Usage:\n\t" << argv[0] << " configuration.yaml" << endl;
    return 1;
  }

  char filename[64], output_folder[64];

  Configuration conf(argv[1]);
  conf.display();
  if (conf.randomise_seed_number < 0) { long _seed = time(NULL); srand(_seed); cout << "Using " << _seed << " to seed a random generator" << endl; }
  if (conf.randomise_seed_number > 0) { long _seed = conf.randomise_seed_number; srand(_seed); cout << "Using " << _seed << " to seed a random generator" << endl; }
  
  int ul_bytes = conf.nMicrolenses * sizeof(Microlens);
  int ray_bytes = conf.nRays * sizeof(Ray);
  int image_bytes = conf.image_height * conf.image_width * sizeof(int);

  Microlens *microlenses = (Microlens*)malloc(ul_bytes);
  Ray *rays = (Ray*)malloc(ray_bytes);
  int *image = (int*)malloc(image_bytes);

  Microlens *ul_buf;
  int *image_buf;
  Ray *ray_buf;

  struct stat info;

  sprintf(output_folder, "./output/%s", conf.configuration_id.c_str());
  if( stat(output_folder , &info ) != 0 ) {
    if (mkdir(output_folder, 0755) != 0 && errno != EEXIST) {
      cerr << "Failed to create output folder " << output_folder << endl;
      return -1;
    }
  }
  
  cout << "Creating microlensing field ... " << flush;
  StartTimer();
  randomiseMicrolenses(microlenses, conf.nMicrolenses, conf.R_field);  
  cudaMalloc(&ul_buf, ul_bytes);
  cudaMemcpy(ul_buf, microlenses, ul_bytes, cudaMemcpyHostToDevice);
  cout << GetElapsedTime() << " s" << endl;
  
  cout << "Defining rays field in ... " << flush;
  StartTimer();
  populateRays(rays, conf.nRays, conf.R_rays, conf.dx_rays);
  cudaMalloc(&image_buf, image_bytes);
  cudaMalloc(&ray_buf, ray_bytes);
  cudaMemcpy(ray_buf, rays, ray_bytes, cudaMemcpyHostToDevice);
  cout << GetElapsedTime() << " s" << endl;

  int nBlocks = (conf.nRays + CUDA_BLOCK_SIZE - 1) / CUDA_BLOCK_SIZE;

  //From https://stackoverflow.com/questions/11888772/when-to-call-cudadevicesynchronize
  //  kernel1<<<X,Y>>>(...); // kernel start execution, CPU continues to next statement
  //  kernel2<<<X,Y>>>(...); // kernel is placed in queue and will start after kernel1 finishes, CPU continues to next statement
  //  cudaMemcpy(...); // CPU blocks until memory is copied, memory copy starts only after kernel2 finishes

  ofstream outf;
  int counter = 0;
  for (float t = 0; t <= conf.t_max; t = t + conf.dt) {
    
    memset(image, 0, image_bytes);
    cudaMemcpy(image_buf, image, image_bytes, cudaMemcpyHostToDevice);

    cout << endl << ">>> Iteration #" << ++counter << ", t=" << t << endl;
    cout << "    Executing ray tracing ... " << flush;
    StartTimer();
    deflectRays<<<nBlocks, CUDA_BLOCK_SIZE>>>(ul_buf, ray_buf, conf, t, image_buf); // compute ray deflections
    //deflectRaysCPU(microlenses, rays, conf, t); // CPU version
    cudaMemcpy(rays, ray_buf, ray_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(image, image_buf, image_bytes, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cout << GetElapsedTime() << " s" << endl;
    
    //sprintf(filename, "%s/rays_y_%.2f.dat", output_folder, t);
    //cout << "  Writing data to " << filename << " ... " << flush;
    //outf.open(filename);
    //for (int i = 0; i <= conf.nRays; i++) {
    //  if (rays[i].d <= conf.image_diagonal_size && rays[i].x1 >= -30 && rays[i].x1 <= 30 && rays[i].x2 >= -30 && rays[i].x2 <= 30) {
    //    outf << rays[i].x1 << " " << rays[i].x2 << endl;
    //  }
    //}
    //outf.close();
    //cout << GetElapsedTime() << " s" << endl;

    sprintf(filename, "%s/image_%.2f.dat", output_folder, t);
    if (conf.save_images) {
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
      cout << GetElapsedTime() << " s" << endl;
    } else {
      cout << "    Skipping writing data to " << filename << " ... " << endl;
    }
  }

  free(microlenses);
  free(rays);
  free(image);
  
  cudaFree(ul_buf);
  cudaFree(ray_buf);
  cudaFree(image_buf);

  return 0;
}
