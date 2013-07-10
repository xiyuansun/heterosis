#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <thrust/reduce.h>

__global__ void sampleTau_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g = GENE, G = a->G;
  
  if(g < G)
    a->tmp1[g] = 1/pow(a->eta[iG(a->mEta, g)], 2);
}

void sampleTau_kernel2(Chain *a){ /* kernel<<<1, 1>>> */
  num_t rate = a->s1 * a->d[a->mD] / 2 + a->bTau;
  num_t shape = a->aTau + a->G * a->d[a->mD] / 2;

  if(shape >= 1 && rate > 0){
    a->tau[a->mTau + 1] = 1/sqrt(rgamma(shape, rate, 0));
  } else {
    a->tau[a->mTau + 1] = a->tau[a->mTau];
  }

  ++a->mTau;
}

void sampleTau(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */
  if(cfg->constTau)
    return;

  __global__ sampleTau_kernel1<<<NBLOCKS, NTHREADS>>>(dev_a);
  
  thrust::device_ptr<num_t> tmp1(host_a->tmp1);  
  num_t s1 = thrust::reduce(tmp1, tmp1 + cfg->G);
  CUDA_CALL(cudaMemcpy(&(dev_a->s1), &s1, sizeof(num_t), cudaMemcpyHostToDevice));
  
  printf("thrust tau s1 = %0.3f\n", s1);
  
  sampleTau_kernel2<<<1, 1>>>(dev_a);
}