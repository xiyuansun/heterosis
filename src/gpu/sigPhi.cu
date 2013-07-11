#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <thrust/reduce.h>

__global__ void sampleSigPhi_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g = IDX, G = a->G;

  if(g < G) 
    a->tmp1[g] = pow(a->phi[iG(a->mPhi, g)] - a->thePhi[a->mThePhi], 2);
}

__global__ void sampleSigPhi_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  num_t rate = a->s1 / 2;
  num_t shape = (a->G - 1) / 2;
  num_t lb = 1/pow(a->sigPhi0, 2);

  if(shape >= 1 && rate > 0){
    a->sigPhi[a->mSigPhi + 1] = 1/sqrt(rgammaDevice(a, 1, shape, rate, lb));
  } else {
    a->sigPhi[a->mSigPhi + 1] = a->sigPhi[a->mSigPhi];
  }
 
  ++a->mSigPhi;
}

void sampleSigPhi(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */
  fprintf(cfg->log, "sigPhi ");

  if(cfg->constSigPhi)
    return;

  sampleSigPhi_kernel1<<<G_GRID, G_BLOCK>>>(dev_a);
  
  thrust::device_ptr<num_t> tmp1(host_a->tmp1);  
  num_t s1 = thrust::reduce(tmp1, tmp1 + cfg->G);
  CUDA_CALL(cudaMemcpy(&(dev_a->s1), &s1, sizeof(num_t), cudaMemcpyHostToDevice));
  
  sampleSigPhi_kernel2<<<1, 1>>>(dev_a);
}