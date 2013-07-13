#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

__global__ void sampleSigC_kernel(Chain *a){ /* kernel <<<1, 1>>> */
  int n, N = a->N;
  num_t rate, shape, lb;

  if(a->constSigC)
    return; 

  rate = 0;
  for(n = 0; n < a->N; ++n) 
    rate += a->c[iN(a->mC, n)] * a->c[iN(a->mC, n)];
  
  shape = (a->N - 1) / 2; 
  rate = rate / 2;
  lb = 1 / pow(a->sigC0, 2); 

  if(shape >= 1 && rate > 0){
    a->sigC[a->mSigC + 1] = 1/sqrt(rgammaDevice(a, 1, shape, rate, lb));
  } else {
    a->sigC[a->mSigC + 1] = a->sigC[a->mSigC];
  }

  ++a->mSigC;
}

__host__ void sampleSigC(Chain *host_a, Chain *dev_a, Config *cfg){ 

  float myTime;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  if(cfg->verbose)
    printf("sigC "); 

  sampleSigC_kernel<<<1, 1>>>(dev_a);
  
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  
  fprintf(cfg->time, "%0.3f ", myTime/MILLISECS); /* elapsed time */
}