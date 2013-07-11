#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <thrust/reduce.h>

__global__ void samplePiAlp_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g = IDX, G = a->G;

  if(g < G){ 
    if(pow(a->alp[iG(a->mAlp, g)], 2) > 1e-6){
      a->tmp1[g] = 1;
    } else {
      a->tmp1[g] = 0;
    }
  }
}

__global__ void samplePiAlp_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  a->piAlp[a->mPiAlp + 1] = rbetaDevice(a, 1, a->G + a->s1 + a->aTau, a->s1 + a->bTau);
  ++a->mPiAlp;
} 

__host__ void samplePiAlp(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */
  fprintf(cfg->log, "piAlp ");

  if(cfg->constPiAlp)
    return;

  samplePiAlp_kernel1<<<G_GRID, G_BLOCK>>>(dev_a);
  
  thrust::device_ptr<num_t> tmp1(host_a->tmp1);  
  num_t s1 = thrust::reduce(tmp1, tmp1 + cfg->G);
  CUDA_CALL(cudaMemcpy(&(dev_a->s1), &s1, sizeof(num_t), cudaMemcpyHostToDevice));
  
  samplePiAlp_kernel2<<<1, 1>>>(dev_a);
}