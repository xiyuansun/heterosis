#include <Chain.h>
#include <constants.h>
#include <constants.h>
#include <cuda.h>
#include <curand_kernel.h>
#include <deviceFunctions.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <thrust/reduce.h>

__global__ void lC_kernel1(Chain *a, int n){ /* kernel <<<G, 1>>> */
  int g = GENE, N = a->N, G = a->G;
  
  if(g < G)
    a->tmp1[g] = exp(a->eps[iNG(a->mEps, n, g)] + mu(a, n, a->phi[iG(a->mPhi, g)], 
                    a->alp[iG(a->mAlp, g)], a->del[iG(a->mDel, g)]));
}

__global__ void lC_kernel2(Chain *a, int n, int newArg){ /* kernel <<<1, 1>>> */
  num_t arg, ret;

  if(newArg){
    arg = a->New[n];
  } else {
    arg = a->Old[n];
  }

  ret = arg * a->G * a->yMeanG[n] - exp(arg) * a->s1 - (arg*arg) / 
        (2 * a->sigC[a->mSigC] * a->sigC[a->mSigC]);

  if(newArg){
    a->lNew[n] = ret;
  } else {
    a->lOld[n] = ret;
  }
}

__host__ void lC(Chain *host_a, Chain *dev_a, Config *cfg, int n, int newArg){ /* host */
  lC_kernel1<<<NBLOCKS, NTHREADS>>>(dev_a, n);
  
  thrust::device_ptr<num_t> tmp1(host_a->tmp1);  
  num_t s1 = thrust::reduce(tmp1, tmp1 + cfg->G);
  CUDA_CALL(cudaMemcpy(&(dev_a->s1), &s1, sizeof(num_t), cudaMemcpyHostToDevice));
  
    printf("thrust c s1 = %0.3f\n", s1);
  
  lC_kernel2<<<1, 1>>>(dev_a, n, newArg);
}

__global__ void sampleC_kernel1(Chain *a){ /* kernel <<<1, N>>> */
  int n = ((blockDim.x * blockIdx.x) + threadIdx.x), N = a->N;
  
  if(n < N){
    a->Old[n] = a->c[iN(a->mC, n)];
    a->New[n] = rnormalDevice(a, n, a->Old[n], a->tuneC[n]);
    
  }
}

__global__ void sampleC_kernel2(Chain *a){ /* kernel <<<1, N>>> */
  int n = ((blockDim.x * blockIdx.x) + threadIdx.x), N = a->N;
  num_t dl, lp, lu;

  if(n < N){ 

    dl = a->lNew[n] - a->lOld[n];
    lp = 0 < dl ? 0 : dl;
    lu = log(runiformDevice(a, 1, 0, 1));
      
    if(lu < lp){ /* accept */
      a->c[iN(a->mC + 1, n)] = a->New[n];
      a->tuneC[n] *= 1.1; /* Increase the proposal variance to avoid  
                                       gettiG stuck in a mode */
                                       
      if(a->mC >= a->burnin)                                 
        ++a->accC[n];
    } else { /* reject */
      a->c[iN(a->mC + 1, n)] = a->Old[n];
      a->tuneC[n] /= 1.1; /* If you're rejectiG too often, decrease the  
                                       proposal variance to sample closer to 
                                       the last accepted value. */
    }
  }
}
 
__global__ void sampleC_kernel3(Chain *a){ /* kernel <<<1, 1>>> */
  ++a->mC;
}

__host__ void sampleC(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */
  int n, N = cfg->N;
  int nthreads = (N < MAXTHREADS ? N : MAXTHREADS);
  int nblocks = ceil(cfg->N / NTHREADS) + 1;
  
  sampleC_kernel1<<<nblocks, nthreads>>>(dev_a);

  for(n = 0; n < cfg->N; ++n){ 
    lC(host_a, dev_a, cfg, n, 1);
    lC(host_a, dev_a, cfg, n, 0);
  }

  sampleC_kernel2<<<nblocks, nthreads>>>(dev_a);
  sampleC_kernel3<<<1, 1>>>(dev_a);
}