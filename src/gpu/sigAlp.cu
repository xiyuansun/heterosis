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

__global__ void sampleSigAlp_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g = IDX, G = a->G;

  if(g < G){
    if(pow(a->alp[iG(a->mAlp, g)], 2) > 1e-6){
      a->tmp1[g] = pow(a->alp[iG(a->mAlp, g)] - a->theAlp[a->mTheAlp], 2);
      a->tmp2[g] = 1;
    } else {
      a->tmp1[g] = 0;
      a->tmp2[g] = 0;
    }
  }
}

__global__ void sampleSigAlp_kernel2(Chain *a){ /* kernel<<<1, 1>>> */
  num_t shape = (a->s2 - 1) / 2;
  num_t rate = a->s1 / 2;
  num_t lb = 1/pow(a->sigAlp0, 2);

  if(shape >= 1 && rate > 0){
    a->sigAlp[a->mSigAlp + 1] = 1/sqrt(rgammaDevice(a, 1, shape, rate, lb));
  } else {
    a->sigAlp[a->mSigAlp + 1] = a->sigAlp[a->mSigAlp]; 
  }

  ++a->mSigAlp; 
}

__host__ void sampleSigAlp(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */
  float myTime;
  cudaError_t err;
  
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  err = cudaEventRecord(start, 0);
if(err != cudaSuccess) {
          printf ("\n\n 1. Error: %s\n\n", cudaGetErrorString(err));
          exit(1);
        }

  fprintf(cfg->log, "sigAlp ");

  if(cfg->constSigAlp)
    return;

  sampleSigAlp_kernel1<<<G_GRID, G_BLOCK>>>(dev_a);
  
  thrust::device_ptr<num_t> tmp1(host_a->tmp1);  
  num_t s1 = thrust::reduce(tmp1, tmp1 + cfg->G);
  CUDA_CALL(cudaMemcpy(&(dev_a->s1), &s1, sizeof(num_t), cudaMemcpyHostToDevice));
  
  thrust::device_ptr<num_t> tmp2(host_a->tmp2);  
  num_t s2 = thrust::reduce(tmp2, tmp2 + cfg->G);
  CUDA_CALL(cudaMemcpy(&(dev_a->s2), &s2, sizeof(num_t), cudaMemcpyHostToDevice));
 
  sampleSigAlp_kernel2<<<1, 1>>>(dev_a); 

  err= cudaEventRecord(stop, 0);
  f(err != cudaSuccess) {
          printf ("\n\n 1. Error: %s\n\n", cudaGetErrorString(err));
          exit(1);
        }
  
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  
  fprintf(cfg->time, "%0.3f ", myTime/MILLISECS); /* elapsed time */
  cudaDeviceSynchronize();
}