#include <common.cuh>
#include <config.cuh>

float distance(float x, float y, float center_x, float center_y) {
  return sqrt(pow(x - center_x, 2) + pow(y - center_y, 2));
}

float distance(float x, float y) {
  return distance(x, y, 0, 0);
}

void randomiseMicrolenses(Microlens *ul, Configuration conf) {
  int n = conf.nMicrolenses;
  float R = conf.R_field;
  float M = conf.M_max;
  float V = conf.V_max;

  for (int i = 0; i < n; i++) {
    float x1 = 2 * R * (rand() / (float)RAND_MAX) - R;
    float x2 = 2 * R * (rand() / (float)RAND_MAX) - R;

    while (distance(x1, x2) > R) {
      x1 = 2 * R * (rand() / (float)RAND_MAX) - R;
      x2 = 2 * R * (rand() / (float)RAND_MAX) - R;
    }
    ul[i] = {.x1 = x1, .x2 = x2, .v1 = 0.0, .v2 = 0.0, .m = 1.0 };
  }

  for (int i = 0; i < n; i++) {
    float v1 = V * (rand() / (float)RAND_MAX) - V;
    float v2 = V * (rand() / (float)RAND_MAX) - V;
    ul[i].v1 = v1;
    ul[i].v2 = v2;
  }

  for (int i = 0; i < n; i++) {
    float m = M * (rand() / (float)RAND_MAX);
    ul[i].m = m;
    //cout << ul[i].x1 << "\t" << ul[i].x2 << "\t" << ul[i].v1 << "\t" << ul[i].v2 << "\t" << ul[i].m << endl;
  }
}

void populateRays(Ray *rays, int nRays, float R_rays, float dx_rays) {
  //for (int i = 0; i < nRays; i++) rays[i] = { .x1 = 0, .x2 = 0 };
  int counter = 0;
  for (float x1 = - R_rays; x1 <= R_rays; x1 += dx_rays) {
    for (float x2 = - R_rays; x2 <= R_rays; x2 += dx_rays) {
      if (distance(x1, x2) <= R_rays && counter < nRays) rays[counter++] = {.x1 = x1, .x2 = x2 };
    }
  }
}

float deg2rad(float deg) {
  return deg * M_PI / 180.0; 
}

void createLC(float *lc_trajectory, const Configuration conf) {
    int counter = 0;
    float y1_cos = cos(deg2rad(conf.lc_angle));
    float y2_sin = sin(deg2rad(conf.lc_angle));
    float y1, y2;
    for (float t = 0.0; t < conf.lc_dist_max; t = t + conf.lc_dist_step) {
          y1 = conf.lc_start_y1 + t * y1_cos;
          y2 = conf.lc_start_y2 + t * y2_sin;
          lc_trajectory[counter + 0 * conf.nLCsteps] = y1; // Y1 coordinate
          lc_trajectory[counter + 1 * conf.nLCsteps] = y2; // Y2 coordinate
          for (int i = 2; i < conf.nLCcolumns; i += 2) {
            lc_trajectory[counter + (i + 0) * conf.nLCsteps] = 0.0; // Amplitude value
            lc_trajectory[counter + (i + 1) * conf.nLCsteps] = 0.0; // Normalization value
          }
        counter++;
    }
}

void resetLC(float *lc_trajectory, const Configuration conf) {
  for (int counter = 0; counter < conf.nLCsteps; counter++) {
    for (int i = 2; i < conf.nLCcolumns; i += 2) {
      lc_trajectory[counter + (i + 0) * conf.nLCsteps] = 0.0; // Amplitude value
      lc_trajectory[counter + (i + 1) * conf.nLCsteps] = 0.0; // Normalization value
    }
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

__device__ float H(float x) {
  return (x > 0) ? 1.0 : 0.0;
}

__global__ void deflectRays(Microlens *uls, Ray *rays, const Configuration c, const float t, int *image, float* lc) {
  int ri = blockDim.x * blockIdx.x + threadIdx.x;
  int j = ri / c.nRaysLine;
  int i = ri - j * c.nRaysLine;
  float ray_x1 = i * c.dx_rays - c.R_rays;
  float ray_x2 = j * c.dx_rays - c.R_rays;
  //if (ri < c.nRays) {
  //  double ray_x1 = rays[ri].x1;
  //  double ray_x2 = rays[ri].x2;

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
    //
    //if (c.save_rays) {
    //    rays[ri].x1 = ray_x1;
    //    rays[ri].x2 = ray_x2;    
    //}
    if (ray_x1 >= c.image_y1_left && ray_x1 <= c.image_y1_right) {
      if (ray_x2 >= c.image_y2_bottom && ray_x2 <= c.image_y2_top) {
        int w = ((ray_x1 - c.image_y1_left) / c.image_pixel_y1_size);
        int h = ((ray_x2 - c.image_y2_bottom) / c.image_pixel_y2_size);
        if (w >= 0 && h >= 0 && w < c.image_width && h < c.image_height) atomicAdd(&image[w * c.image_height + h], 1);
      }
    }
//    int w = ((ray_x1 - c.image_y1_left) / c.image_pixel_y1_size);
//    int h = ((ray_x2 - c.image_y2_bottom) / c.image_pixel_y2_size);
//    if (w >= 0 && h >= 0 && w < c.image_width && h < c.image_height) atomicAdd(&image[w * c.image_height + h], 1);
  }
}

__global__ void calculateLCs(const Configuration c, int *image, float *lc) {
  int w = blockIdx.x * blockDim.x + threadIdx.x;
  int h = blockIdx.y * blockDim.y + threadIdx.y;
  int pix = image[w * c.image_height + h];

  if ((pix > 0) && (w * c.image_height + h < c.image_height * c.image_width)) {
    float r_x1 = w * c.image_pixel_y1_size + c.image_y1_left;
    float r_x2 = h * c.image_pixel_y2_size + c.image_y2_bottom;
  
    float factorex_gs, factorex_ld, factorex_pl, factorex_ad, factorex_el, factorex_el_p;
    for (int i = 0; i < c.nLCsteps - 1; i++) {
      float lc_y1 = lc[i + 0 * c.nLCsteps];
      float lc_y2 = lc[i + 1 * c.nLCsteps];
      float d = dst(lc_y1 - r_x1, lc_y2 - r_x2);
      float d2 = d * d;
 
      int index = 2;
      for (float source_size = c.source_size[0]; source_size < c.source_size[1]; source_size+=c.source_size[2]) {
        if (d < 10 * source_size) {

          float R_gs = source_size;
          float R2_gs = R_gs * R_gs;
      
          float R_1_2_ld = source_size * sqrt(log(2.0));
          float R_ld = R_1_2_ld / sqrt(1.0 - pow(0.5, 2)/(c.p_ld + 2));
          float R2_ld = R_ld * R_ld;
    
          float R_1_2_pl = source_size * sqrt(log(2.0));
          float R_pl = R_1_2_pl / sqrt((pow(2.0, 1.0/(c.p_pl - 1)) - 1.0)/log(2.0));
          float R2_pl = R_pl * R_pl;
      
          float R_1_2_ad = source_size * sqrt(log(2.0));
          float R_ad = R_1_2_ad/4.0;
          //float R2_ad = R_ad * R_ad; // unused
    
          factorex_gs = expf(- d2 / R2_gs);
          factorex_ld = ((c.p_ld + 1)/(M_PI * R2_ld)) * H(1 - d2/R2_ld) * pow(1 - d2/R2_ld, c.p_ld);
          factorex_pl = ((c.p_pl - 1)/(M_PI * R2_pl)) * (1/pow(1 + d2/R2_pl, c.p_pl));
          if (d > R_ad) {
            factorex_ad = (3 * R_ad  / (2 * M_PI * pow(d, 3))) * (1 - sqrt(R_ad/d));
            atomicAdd(&lc[i + (index + 0) * c.nLCsteps], factorex_ad); // Normalization
            atomicAdd(&lc[i + (index + 1) * c.nLCsteps], pix * factorex_ad); // Amplitude value, non-normalized  
          }
          atomicAdd(&lc[i + (index + 2) * c.nLCsteps], factorex_gs); // Normalization
          atomicAdd(&lc[i + (index + 3) * c.nLCsteps], pix * factorex_gs); // Amplitude value, non-normalized  
  
          atomicAdd(&lc[i + (index + 4) * c.nLCsteps], factorex_ld); // Normalization
          atomicAdd(&lc[i + (index + 5) * c.nLCsteps], pix * factorex_ld); // Amplitude value, non-normalized  
  
          atomicAdd(&lc[i + (index + 6) * c.nLCsteps], factorex_pl); // Normalization
          atomicAdd(&lc[i + (index + 7) * c.nLCsteps], pix * factorex_pl); // Amplitude value, non-normalized  
          
          int e_index = index + 8;
          for (float eccentricity = c.eccentricity[0]; eccentricity < c.eccentricity[1]; eccentricity+=c.eccentricity[2]) {
            float e2_el = eccentricity * eccentricity;
            float a_el = source_size / (1 - e2_el);
            float b_el = source_size * (1 - e2_el);
        
            //float a2_el = a_el * a_el; // unused
            //float b2_el = b_el * b_el; // unused

            float d_el = dst((lc_y1 - r_x1) / a_el, (lc_y2 - r_x2) / b_el);
            float d2_el = d_el * d_el;
    
            float d_el_p = dst((lc_y1 - r_x1) / b_el, (lc_y2 - r_x2) / a_el);
            float d2_el_p = d_el_p * d_el_p;

            factorex_el = expf( - d2_el);
            factorex_el_p = expf( - d2_el_p);
  
            atomicAdd(&lc[i + (e_index + 0) * c.nLCsteps], factorex_el); // Normalization
            atomicAdd(&lc[i + (e_index + 1) * c.nLCsteps], pix * factorex_el); // Amplitude value, non-normalized  

            atomicAdd(&lc[i + (e_index + 2) * c.nLCsteps], factorex_el_p); // Normalization
            atomicAdd(&lc[i + (e_index + 3) * c.nLCsteps], pix * factorex_el_p); // Amplitude value, non-normalized  

            lc[(c.nLCsteps - 1) + (e_index + 0) * c.nLCsteps] = eccentricity;
            lc[(c.nLCsteps - 1) + (e_index + 1) * c.nLCsteps] = 5;

            lc[(c.nLCsteps - 1) + (e_index + 2) * c.nLCsteps] = eccentricity;
            lc[(c.nLCsteps - 1) + (e_index + 3) * c.nLCsteps] = 6;

            e_index += 4;
          }
          lc[(c.nLCsteps - 1) + (index + 0) * c.nLCsteps] = R_gs;
          lc[(c.nLCsteps - 1) + (index + 1) * c.nLCsteps] = 1;
  
          lc[(c.nLCsteps - 1) + (index + 2) * c.nLCsteps] = R_gs;
          lc[(c.nLCsteps - 1) + (index + 3) * c.nLCsteps] = 2;
  
          lc[(c.nLCsteps - 1) + (index + 4) * c.nLCsteps] = R_gs;
          lc[(c.nLCsteps - 1) + (index + 5) * c.nLCsteps] = 3;
  
          lc[(c.nLCsteps - 1) + (index + 6) * c.nLCsteps] = R_gs;
          lc[(c.nLCsteps - 1) + (index + 7) * c.nLCsteps] = 4;
        }
        index += 8 + 4 * c.nCountEccentricities;
      }
      /*
      if (d < 10 * c.R_gs) {
        float d_el = dst((lc_y1 - r_x1) / c.a_el, (lc_y2 - r_x2) / c.b_el);
        float d2_el = d_el * d_el;

        float d_el_p = dst((lc_y1 - r_x1) / c.b_el, (lc_y2 - r_x2) / c.a_el);
        float d2_el_p = d_el_p * d_el_p;
          
        factorex_gs = expf(- d2 / c.R2_gs);
        factorex_ld = ((c.p_ld + 1)/(M_PI * c.R2_ld)) * H(1 - d2/c.R2_ld) * pow(1 - d2/c.R2_ld, c.p_ld);
        factorex_pl = ((c.p_pl - 1)/(M_PI * c.R2_pl)) * (1/pow(1 + d2/c.R2_pl, c.p_pl));
        factorex_el = expf( - d2_el);
        factorex_el_p = expf( - d2_el_p);
        if (d > c.R_ad) {
          factorex_ad = (3 * c.R_ad  / (2 * M_PI * pow(d, 3))) * (1 - sqrt(c.R_ad/d));
          atomicAdd(&lc[i + 2 * c.nLCsteps], factorex_ad); // Normalization
          atomicAdd(&lc[i + 3 * c.nLCsteps], pix * factorex_ad); // Amplitude value, non-normalized  
        }
        atomicAdd(&lc[i + 4 * c.nLCsteps], factorex_gs); // Normalization
        atomicAdd(&lc[i + 5 * c.nLCsteps], pix * factorex_gs); // Amplitude value, non-normalized  

        atomicAdd(&lc[i + 6 * c.nLCsteps], factorex_ld); // Normalization
        atomicAdd(&lc[i + 7 * c.nLCsteps], pix * factorex_ld); // Amplitude value, non-normalized  

        atomicAdd(&lc[i + 8 * c.nLCsteps], factorex_pl); // Normalization
        atomicAdd(&lc[i + 9 * c.nLCsteps], pix * factorex_pl); // Amplitude value, non-normalized  

        atomicAdd(&lc[i + 10 * c.nLCsteps], factorex_el); // Normalization
        atomicAdd(&lc[i + 11 * c.nLCsteps], pix * factorex_el); // Amplitude value, non-normalized  

        atomicAdd(&lc[i + 12 * c.nLCsteps], factorex_el_p); // Normalization
        atomicAdd(&lc[i + 13 * c.nLCsteps], pix * factorex_el_p); // Amplitude value, non-normalized  
      }
      */
    }   
  }
}
