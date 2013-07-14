#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

__global__ void sampleSigC_kernel(Chain *a){ /* kernel <<<1, 1>>> */
  int n;
  num_t shape, rate = 0, lb;
  
  for(n = 0; n < a->N; ++n) 
    rate += a->c[n] * a->c[n];
  
  shape = (a->N - 1) / 2; 
  rate = rate / 2;
  lb = 1 / pow((float) a->sigC0, 2); 

  if(shape >= 1 && rate > 0){
    a->sigC = 1/sqrt(rgammaDevice(a, 1, shape, rate, lb));
  } 
}

__host__ void sampleSigC(Chain *host_a, Chain *dev_a, Config *cfg){ 

  num_t myTime;
  cudaEvent_t start, stop;
  
  if(cfg->constSigC)
    return; 
  
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

  cfg->timeSigC = myTime / MILLISECS;
}