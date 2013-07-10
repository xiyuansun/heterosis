#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <thrust/reduce.h>

__global__ void lD_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g = GENE, G = a->G;

  if(g < G){ 
    a->tmp1[g] = 2 * log(a->eta[iG(a->mEta, g)]);
    a->tmp2[g] = 1/(a->eta[iG(a->mEta, g)] * a->eta[iG(a->mEta, g)]);
  }
}

__global__ void lD_kernel2(Chain *a, int newArg){ /* kernel <<<1, 1>>> */
  num_t arg, ret, tmp;
 
  if(newArg){
    arg = a->New[0];
  } else{
    arg = a->Old[0];
  }

  tmp = arg * a->tau[a->mTau] * a->tau[a->mTau] / 2;
  ret = -a->G * lgamma(arg/2) + (a->G * arg / 2) * log(tmp);
  ret -= (arg/2 + 1) * a->s1 - tmp * a->s2;

  if(newArg){
    a->lNew[0] = ret;
  } else{
    a->lOld[0] = ret;
  }
}

__host__ void lD(Chain *host_a, Chain *dev_a, Config *cfg, int newArg){ /* host */
  
  if(newArg){
    if(host_a->New[0] <= 0 || host_a->New[0] > host_a->d0)
     host_a->lNew[0] = NUM_TMIN;
  } else {
    if(host_a->Old[0] <= 0 || host_a->Old[0] > host_a->d0)
      host_a->lOld[0] = NUM_TMIN; 
  }

  lD_kernel1<<<NBLOCKS, NTHREADS>>>(dev_a);
  
  thrust::device_ptr<num_t> tmp1(host_a->tmp1);  
  num_t s1 = thrust::reduce(tmp1, tmp1 + cfg->G);
  CUDA_CALL(cudaMemcpy(&(dev_a->s1), &s1, sizeof(num_t), cudaMemcpyHostToDevice));
  
  thrust::device_ptr<num_t> tmp2(host_a->tmp2);  
  num_t s2 = thrust::reduce(tmp2, tmp2 + cfg->G);
  CUDA_CALL(cudaMemcpy(&(dev_a->s2), &s2, sizeof(num_t), cudaMemcpyHostToDevice));
  
  lD_kernel2<<<1, 1>>>(dev_a, newArg);
}

__global__ void sampleD_kernel1(Chain *a){ /* kernel <<<1, 1>>> */
  a->Old[0] = a->d[a->mD];
  
  do {
    a->New[0] = rnormal(a->Old[0], a->tuneD);
  } while(a->New[0] < 0);
}

__global__ void sampleD_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  num_t dl = a->lNew[0] - a->lOld[0];
  num_t lp = 0 < dl ? 0 : dl;
  num_t lu = log(runiform(0, 1));

  if(lu < lp){ /* accept */
    a->d[a->mD + 1] = a->New[0];
    a->tuneD *= 1.1; /* Increase the proposal variance to avoid getting 
                                  stuck in a mode */
    
    if(a->mD >= a->burnin) 
      ++a->accD;
  } else { /* reject */
    a->d[a->mD + 1] = a->Old[0];
    a->tuneD /= 1.1; /* If you're rejecting too often, decrease the proposal 
                                  variance to sample closer to the last accepted value. */
  }

  ++a->mD;
}

__host__ void sampleD(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */
  if(cfg->constD)
    return;
   
  sampleD_kernel1<<<1, 1>>>(dev_a);

  lD(host_a, dev_a, 1);
  lD(host_a, dev_a, 0);

  sampleD_kernel2<<<1, 1>>>(dev_a);
}
