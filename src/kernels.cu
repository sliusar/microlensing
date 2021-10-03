#include <common.cuh>
#include <config.cuh>

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
    for (int i = 0; i < nRays; i++) rays[i] = { .x1 = 0, .x2 = 0 };
    //int counter = 0;
    //for (float x1 = - R_rays; x1 <= R_rays; x1 += dx_rays) {
    //  for (float x2 = - R_rays; x2 <= R_rays; x2 += dx_rays) {
    //    if (distance(x1, x2) <= R_rays && counter < nRays) rays[counter++] = {.x1 = x1, .x2 = x2 };
    //  }
    //}
}

void createTrajectory(float *lc_trajectory, const Configuration conf) {
    int counter = 0;
    float l = sqrt(pow(conf.lc_end_y1 - conf.lc_start_y1, 2) + pow(conf.lc_end_y2 - conf.lc_start_y2, 2));
    float step_y1 = conf.lc_step * (conf.lc_end_y1 - conf.lc_start_y1)/l;
    float step_y2 = conf.lc_step * (conf.lc_end_y2 - conf.lc_start_y2)/l;
    for (float y1 = conf.lc_start_y1, y2 = conf.lc_start_y2; (y1 < conf.lc_end_y1 && y2 < conf.lc_end_y2); y1+=step_y1, y2+=step_y2) {
        //i_y1 = round((y1 - conf.image_y1_left) / conf.image_pixel_y1_size);
        //i_y2 = round((y2 - conf.image_y2_bottom) / conf.image_pixel_y2_size);
        lc_trajectory[counter + 0 * conf.nLCsteps] = y1; // Y1 coordinate
        lc_trajectory[counter + 1 * conf.nLCsteps] = y2; // Y2 coordinate
        lc_trajectory[counter + 2 * conf.nLCsteps] = 0.0; // Gauss amplitude normalization value
        lc_trajectory[counter + 3 * conf.nLCsteps] = 0.0; // Gauss amplitude value
        counter++;
    }
}

void printLC(float *lc, int c) {
  cout << endl << "Y1\tY2\tNorm\tVal" << endl;
  for (int i = 0; i < c; i++) {
    cout << lc[i + 0 * c] << "\t" << lc[i + 1 * c] << "\t" << lc[i + 2 * c] << "\t" << lc[i + 3 * c] << endl;
  }
}

__device__ float dst2_inv(float x, float y) {
  return rhypotf(x, y) * rhypotf(x, y);
}

__device__ float dst(float x, float y) {
  return hypotf(x, y);
}

__device__ float dst(float x1, float y1, float x2, float y2) {
  return hypotf(x1 - x2, y1 - y2);
}

__device__ float H(float a, float x) {
  return (x >= a) ? 1.0 : 0.0;
}

__global__ void deflectRays(Microlens *uls, Ray *rays, const Configuration c, const float t, int *image, float* lc) {
  int ri = blockDim.x * blockIdx.x + threadIdx.x;
  int j = ri / c.nRays_line;
  int i = ri - j * c.nRays_line;
  float ray_x1 = i * c.dx_rays - c.R_rays;
  float ray_x2 = j * c.dx_rays - c.R_rays;

  if (dst(ray_x1, ray_x2) <= c.R_rays) {
    float sum_x1 = 0.0;
    float sum_x2 = 0.0;
    for (int i = 0; i < c.nMicrolenses; i++) {
      float m_x1 = ray_x1 - uls[i].x1 - (uls[i].v1 * t);
      float m_x2 = ray_x2 - uls[i].x2 - (uls[i].v2 * t);
      float ri = uls[i].m * dst2_inv(m_x1, m_x2);
      sum_x1 += m_x1 * ri;
      sum_x2 += m_x2 * ri;
    }
    ray_x1 = (1 - c.gamma) * ray_x1 - c.sigma_c * ray_x1 - sum_x1;
    ray_x2 = (1 + c.gamma) * ray_x2 - c.sigma_c * ray_x2 - sum_x2;

    if (c.output_rays) {
        rays[ri].x1 = ray_x1;
        rays[ri].x2 = ray_x2;    
    }

    int w = lrintf((ray_x1 - c.image_y1_left) / c.image_pixel_y1_size);
    int h = lrintf((ray_x2 - c.image_y2_bottom) / c.image_pixel_y2_size);
    if (w >= 0 && w < c.image_width && h >= 0 && h < c.image_height) atomicAdd(&image[w * c.image_height + h], 1);


    for (int i = 0; i < c.nLCsteps; i++) {
        float lc_y1 = lc[i + 0 * c.nLCsteps];
        float lc_y2 = lc[i + 1 * c.nLCsteps];
        float d = dst(lc_y1 - ray_x1, lc_y2 - ray_x2);
        float sigma = 0.1;
        float sigma2 = 0.01;
        if (d < 4 * sigma) {
            float val = expf(- d * d / sigma2);
            //atomicAdd(&lc[i + 2 * c.nLCsteps], 1.0); // Normalization
            //atomicAdd(&lc[i + 3 * c.nLCsteps], val); // Amplitude value, non-normalized
        }
    }
  }
}

__global__ void calculateLCs(const Configuration c, int *image, float *lc) {
  int w = blockIdx.x * blockDim.x + threadIdx.x;
  int h = blockIdx.y * blockDim.y + threadIdx.y;


  float r_x1 = w * c.image_pixel_y1_size + c.image_y1_left;
  float r_x2 = h * c.image_pixel_y2_size + c.image_y2_bottom;

  for (int i = 0; i < c.nLCsteps; i++) {
    float lc_y1 = lc[i + 0 * c.nLCsteps];
    float lc_y2 = lc[i + 1 * c.nLCsteps];
    float d = dst(lc_y1 - r_x1, lc_y2 - r_x2);
    float sigma = 0.1;
    float sigma2 = 0.01;
    if (d < 4 * sigma) {
        float v = expf(- d * d / sigma2);
        float val = image[w * c.image_height + h] * v;
        atomicAdd(&lc[i + 2 * c.nLCsteps], v); // Normalization
        atomicAdd(&lc[i + 3 * c.nLCsteps], val); // Amplitude value, non-normalized
    }
    //factorex = exp(-dist2 / sigsq2);
    //dum = pix[ix + i1 + (iy + i2) * pix_dim1];
    //value += dum * factorex;
    //normfac += factorex;
  }

}