#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <thrust/reduce.h>

__global__ void sampleSigAlp_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g = IDX;

  if(g < a->G){
    if(pow((float) a->alp[g], 2) > 1e-6){
      a->tmp1[g] = pow((float) (a->alp[g] - a->theAlp), 2);
      a->tmp2[g] = 1;
    } else {
      a->tmp1[g] = 0;
      a->tmp2[g] = 0;
    }
  }
}

__global__ void sampleSigAlp_kernel2(Chain *a){ /* parallel pairwise sum in Thrust */
  num_t shape = (a->s2 - 1) / 2;
  num_t rate = a->s1 / 2;
  num_t lb = 1/pow((float) a->sigAlp0, 2);

  if(shape > 0 && rate > 0){
    a->sigAlp = 1/sqrt(rgammaDevice(a, 1, shape, rate, lb));
  } 
}

__host__ void sampleSigAlp(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */

  float myTime;
  cudaEvent_t start, stop;
  
  if(cfg->constSigAlp)
    return;
  
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  
  if(cfg->verbose)
    printf("sigAlp ");

  if(!cfg->alpPrior){
	sampleSigAlp_kernel1<<<G_GRID, G_BLOCK>>>(dev_a);
  
	thrust::device_ptr<num_t> tmp1(host_a->tmp1);  
	num_t s1 = thrust::reduce(tmp1, tmp1 + cfg->G);
	CUDA_CALL(cudaMemcpy(&(dev_a->s1), &s1, sizeof(num_t), cudaMemcpyHostToDevice));
  
	thrust::device_ptr<num_t> tmp2(host_a->tmp2);  
	num_t s2 = thrust::reduce(tmp2, tmp2 + cfg->G);
	CUDA_CALL(cudaMemcpy(&(dev_a->s2), &s2, sizeof(num_t), cudaMemcpyHostToDevice));
 
	sampleSigAlp_kernel2<<<1, 1>>>(dev_a); 
  }
  
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  cfg->timeSigAlp = myTime / MILLISECS;
}
