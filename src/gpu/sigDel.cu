#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <thrust/reduce.h>

__global__ void sampleSigDel_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g = IDX;

  if(g < a->G){ 
    if(pow((float) a->del[g], 2) > 1e-6){
      a->tmp1[g] = pow((float) (a->del[g] - a->theDel), 2);
      a->tmp2[g] = 1;
    } else {
      a->tmp1[g] = 0;
      a->tmp2[g] = 0;
    }
  }
}

__global__ void sampleSigDel_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  num_t shape = (a->s2 - 1) / 2;
  num_t rate = a->s1 / 2;
  num_t lb = 1/pow((float) a->sigDel0, 2);

  if(shape > 0 && rate > 0){
    a->sigDel = 1/sqrt(rgammaDevice(a, 1, shape, rate, lb));
  } else {
    a->sigDel = a->sigDel;
  }
}

__host__ void sampleSigDel(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */

  float myTime;
  cudaEvent_t start, stop;
  
  if(cfg->constSigDel || !cfg->heterosis)
    return;
  
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  if(cfg->verbose)
    printf("sigDel ");

  if(!cfg->delPrior){
	sampleSigDel_kernel1<<<G_GRID, G_BLOCK>>>(dev_a);
  
	thrust::device_ptr<num_t> tmp1(host_a->tmp1);  
	num_t s1 = thrust::reduce(tmp1, tmp1 + cfg->G);
	CUDA_CALL(cudaMemcpy(&(dev_a->s1), &s1, sizeof(num_t), cudaMemcpyHostToDevice));
  
	thrust::device_ptr<num_t> tmp2(host_a->tmp2);  
	num_t s2 = thrust::reduce(tmp2, tmp2 + cfg->G);
	CUDA_CALL(cudaMemcpy(&(dev_a->s2), &s2, sizeof(num_t), cudaMemcpyHostToDevice));
 
	sampleSigDel_kernel2<<<1, 1>>>(dev_a);
  }

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  cfg->timeSigDel = myTime / MILLISECS;
}
