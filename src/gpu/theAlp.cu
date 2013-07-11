#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <thrust/reduce.h>

__global__ void sampleTheAlp_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g = IDX, G = a->G;

  if(g < G){
    if(pow(a->alp[iG(a->mAlp, g)], 2) > 1e-6){
      a->tmp1[g] = 1;
      a->tmp2[g] = a->alp[iG(a->mAlp, g)];
    } else {
      a->tmp1[g] = 0;
      a->tmp2[g] = 0;
    }
  }
} 

__global__ void sampleTheAlp_kernel2(Chain *a){ /* kernel <<<1, 1>>> */

  num_t gs = a->gamAlp * a->gamAlp;
  num_t ss = a->sigAlp[a->mSigAlp] * a->sigAlp[a->mSigAlp];
  num_t den = a->s1 * gs + ss;

  num_t m = gs * a->s2 / den;
  num_t s = sqrt(gs * ss / den);

  a->theAlp[a->mTheAlp + 1] = rnormalDevice(a, 1, m, s);
  ++a->mTheAlp;
}

__host__ void sampleTheAlp(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */
  fprintf(cfg->log, "theAlp ");

  if(cfg->constTheAlp)
    return;

  sampleTheAlp_kernel1<<<G_GRID, G_BLOCK>>>(dev_a);
  
  thrust::device_ptr<num_t> tmp1(host_a->tmp1);  
  num_t s1 = thrust::reduce(tmp1, tmp1 + cfg->G);
  CUDA_CALL(cudaMemcpy(&(dev_a->s1), &s1, sizeof(num_t), cudaMemcpyHostToDevice));
  
  thrust::device_ptr<num_t> tmp2(host_a->tmp2);  
  num_t s2 = thrust::reduce(tmp2, tmp2 + cfg->G);
  CUDA_CALL(cudaMemcpy(&(dev_a->s2), &s2, sizeof(num_t), cudaMemcpyHostToDevice));
  
  sampleTheAlp_kernel2<<<G_GRID, G_BLOCK>>>(dev_a);
}
