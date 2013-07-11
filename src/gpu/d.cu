#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <thrust/reduce.h>

__global__ void lD_kernel1(Chain *a, int newArg){ /* kernel <<<1, 1>>> */
  if(newArg){
    if(a->New[0] <= 0 || a->New[0] > a->d0)
     a->lNew[0] = NUM_TMIN;
  } else {
    if(a->Old[0] <= 0 || a->Old[0] > a->d0)
      a->lOld[0] = NUM_TMIN; 
  }
}

__global__ void lD_kernel2(Chain *a){ /* kernel <<<G, 1>>> */
  int g = IDX, G = a->G;

  if(g < G){ 
    a->tmp1[g] = 2 * log(a->eta[iG(a->mEta, g)]);
    a->tmp2[g] = 1/(a->eta[iG(a->mEta, g)] * a->eta[iG(a->mEta, g)]);
  }
}

__global__ void lD_kernel3(Chain *a, int newArg){ /* kernel <<<1, 1>>> */
  num_t arg, ret, tmp;
 
  if(newArg){ 
    arg = a->New[0];
  } else{
    arg = a->Old[0];
  }

  tmp = arg * a->tau[a->mTau] * a->tau[a->mTau] / 2;
  ret = -a->G * lgamma(arg/2) + (a->G * arg / 2) * log(tmp);
  ret -= (arg/2 + 1) * a->s1 - tmp * a->s2;

  if(newArg){
    a->lNew[0] = ret;
  } else{
    a->lOld[0] = ret;
  }
}

__host__ void lD(Chain *host_a, Chain *dev_a, Config *cfg, int newArg){ /* host */
  
  lD_kernel1<<<1, 1>>>(dev_a, newArg);
  lD_kernel2<<<G_GRID, G_BLOCK>>>(dev_a);
  
  thrust::device_ptr<num_t> tmp1(host_a->tmp1);  
  num_t s1 = thrust::reduce(tmp1, tmp1 + cfg->G);
  CUDA_CALL(cudaMemcpy(&(dev_a->s1), &s1, sizeof(num_t), cudaMemcpyHostToDevice));
  
  thrust::device_ptr<num_t> tmp2(host_a->tmp2);  
  num_t s2 = thrust::reduce(tmp2, tmp2 + cfg->G);
  CUDA_CALL(cudaMemcpy(&(dev_a->s2), &s2, sizeof(num_t), cudaMemcpyHostToDevice));
  
  lD_kernel3<<<1, 1>>>(dev_a, newArg);
}

__global__ void sampleD_kernel1(Chain *a){ /* kernel <<<1, 1>>> */
  a->Old[0] = a->d[a->mD];
  
  do {
    a->New[0] = rnormalDevice(a, 1, a->Old[0], a->tuneD);
  } while(a->New[0] < 0);
}

__global__ void sampleD_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  num_t dl = a->lNew[0] - a->lOld[0];
  num_t lp = 0 < dl ? 0 : dl;
  num_t lu = log(runiformDevice(a, 1, 0, 1));

  if(lu < lp){ /* accept */
    a->d[a->mD + 1] = a->New[0];
    a->tuneD *= 1.1; /* Increase the proposal variance to avoid getting 
                                  stuck in a mode */
    if(a->mD >= a->burnin) 
      ++a->accD;
  } else { /* reject */
    a->d[a->mD + 1] = a->Old[0];
    a->tuneD /= 1.1; /* If you're rejecting too often, decrease the proposal 
                                  variance to sample closer to the last accepted value. */
  }

  ++a->mD;
}

__host__ void sampleD(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */

  float myTime;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  fprintf(cfg->log, "d ");

  if(cfg->constD)
    return;
   
  sampleD_kernel1<<<1, 1>>>(dev_a);

  lD(host_a, dev_a, cfg, 1);
  lD(host_a, dev_a, cfg, 0);

  sampleD_kernel2<<<1, 1>>>(dev_a);

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  
  fprintf(cfg->time, "%0.3f ", myTime/MILLISECS); /* elapsed time */
}
