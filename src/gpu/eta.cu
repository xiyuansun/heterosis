#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

__global__ void sampleEta_kernel1(Chain *a){ /* kernel <<<1, 1>>> */
  a->s1 = (a->N + a->d) / 2; 
}

__global__ void sampleEta_kernel2(Chain *a){ /* kernel <<<G, 1>>> */
  int n, g = IDX, G = a->G;
  num_t shape = a->s1, rate;

  if(g < G){

    rate = 0;
    for(n = 0; n < a->N; ++n) 
      rate += pow(a->eps[iG(n, g)], 2);
  
    rate = (rate + a->d * a->tau * a->tau) / 2; 

    if(shape >= 1 && rate > 0){
      a->eta[g] = 1/sqrt(rgammaDevice(a, g, shape, rate, 0));
    } else {
      a->eta[g] = a->eta[g];
    }
  }
}

__host__ void sampleEta(Chain *host_a, Chain *dev_a, Config *cfg){

  float myTime;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  if(cfg->verbose)
    printf("eta ");

  sampleEta_kernel1<<<1, 1>>>(dev_a);
  sampleEta_kernel2<<<G_GRID, G_BLOCK>>>(dev_a);
  
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  
  cfg->timeEta = myTime;
}