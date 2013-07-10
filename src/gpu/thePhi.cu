#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <thrust/reduce.h>

__global__ void sampleThePhi_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g = GENE, G = a->G;
  
  if(g < G)
    a->tmp1[g] = a->phi[iG(a->mPhi, g)];
}

__global__ void sampleThePhi_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  num_t gs = a->gamPhi * a->gamPhi;
  num_t ss = a->sigPhi[a->mSigPhi] * a->sigPhi[a->mSigPhi];
  num_t den = (a->G * gs + ss);

  num_t m = gs * a->s1 / den;
  num_t s = gs * ss / den;

  a->thePhi[a->mThePhi + 1] = rnormalDevice(a, 1, m, s);
  ++a->mThePhi;
}

__host__ void sampleThePhi(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */
  if(cfg->constThePhi)
    return;

  sampleThePhi_kernel1<<<NBLOCKS, NTHREADS>>>(dev_a);
  
  thrust::device_ptr<num_t> tmp1(host_a->tmp1);  
  num_t s1 = thrust::reduce(tmp1, tmp1 + cfg->G);
  CUDA_CALL(cudaMemcpy(&(dev_a->s1), &s1, sizeof(num_t), cudaMemcpyHostToDevice));
  
  
  sampleThePhi_kernel2<<<1, 1>>>(dev_a);
}