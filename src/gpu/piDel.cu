#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <thrust/reduce.h>

__global__ void samplePiDel_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g = IDX;

  if(g < a->G){ 
    if(pow((float) a->del[g], 2) > 1e-6){
      a->tmp1[g] = 1; 
    } else {
      a->tmp1[g] = 0;
    }
  } 
}

__global__ void samplePiDel_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  a->piDel = rbetaDevice(a, 1, a->G + a->s1 + a->aTau, a->s1 + a->bTau);
}

__host__ void samplePiDel(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */

  float myTime;
  cudaEvent_t start, stop;
  
  if(cfg->constPiDel || !cfg->heterosis)
    return;
  
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  if(cfg->verbose)
    printf("piDel ");

  samplePiDel_kernel1<<<G_GRID, G_BLOCK>>>(dev_a);
  
  thrust::device_ptr<num_t> tmp1(host_a->tmp1);  
  num_t s1 = thrust::reduce(tmp1, tmp1 + cfg->G);
  CUDA_CALL(cudaMemcpy(&(dev_a->s1), &s1, sizeof(num_t), cudaMemcpyHostToDevice));
  
  samplePiDel_kernel2<<<1, 1>>>(dev_a);

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  cfg->timePiDel = myTime / MILLISECS;
}
