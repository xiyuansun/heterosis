#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <thrust/reduce.h>

__global__ void sampleThePhi_kernel1(Chain *a){ /* pairwise sum in Thrust */
  int g = IDX;
  
  a->s1 = 0; 
  if(g < a->G)
    a->s1 += a->phi[g];
}

__global__ void sampleThePhi_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  num_t gs = a->gamPhi * a->gamPhi;
  num_t ss = a->sigPhi * a->sigPhi;
  num_t den = (a->G * gs + ss);

  num_t m = gs * a->s1 / den;
  num_t s = gs * ss / den;

  a->thePhi = rnormalDevice(a, 1, m, s);
}

__host__ void sampleThePhi(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */

  num_t myTime;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  if(cfg->verbose)
    printf("thePhi ");

  if(cfg->constThePhi)
    return;

  sampleThePhi_kernel1<<<G_GRID, G_BLOCK>>>(dev_a);
  
  thrust::device_ptr<num_t> tmp1(host_a->tmp1);  
  num_t s1 = thrust::reduce(tmp1, tmp1 + cfg->G);
  CUDA_CALL(cudaMemcpy(&(dev_a->s1), &s1, sizeof(num_t), cudaMemcpyHostToDevice));
  
  sampleThePhi_kernel2<<<1, 1>>>(dev_a);

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop); 
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  
  cfg->timeThePhi = myTime;
}