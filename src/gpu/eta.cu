#include <Chain.h>
#include <constants.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <thrust/reduce.h>

__global__ void sampleEta_kernel1(Chain *a){ /* kernel <<<1, 1>>> */
  a->s1 = (a->N + a->d[a->mD]) / 2; 
}

__global__ void sampleEta_kernel2(Chain *a){ /* kernel <<<G, 1>>> */
  int n, g = IDX, N = a->N, G = a->G;
  num_t shape = a->s1, rate;

  if(g < G){

    rate = 0;
    for(n = 0; n < a->N; ++n) 
      rate += pow(a->eps[iNG(a->mEps, n, g)], 2);
  
    rate = (rate + a->d[a->mD] * a->tau[a->mTau] * a->tau[a->mTau]) / 2; 

    if(shape >= 1 && rate > 0){
      a->eta[iG(a->mEta + 1, g)] = 1/sqrt(rgammaDevice(a, 1, shape, rate, 0));
    } else {
      a->eta[iG(a->mEta + 1, g)] = a->eta[iG(a->mEta, g)];
    }
  }
}

__global__ void sampleEta_kernel3(Chain *a){ /* kernel <<<1, 1>>> */
  ++a->mEta;
}

void sampleEta(Chain *host_a, Chain *dev_a, Config *cfg){
  float myTime;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  fprintf(cfg->log, "eta ");

  sampleEta_kernel1<<<1, 1>>>(dev_a);
  sampleEta_kernel2<<<G_GRID, G_BLOCK>>>(dev_a);
  sampleEta_kernel3<<<1, 1>>>(dev_a);
  
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  
  fprintf(cfg->time, "%0.3f ", myTime/MILLISECS); /* elapsed time */
  cudaDeviceSynchronize();
}