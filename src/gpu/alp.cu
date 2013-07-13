#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

__device__ num_t lAlp(Chain *a, int g, num_t arg){ /* device */
  
  int n, N = a->N, G = a->G;
  num_t s = 0, tmp;
   
  for(n = 0; n < a->N; ++n){
    if(a->grp[n] != 2){
      tmp = mu(a, n, a->phi[iG(a->mPhi, g)], arg, a->del[iG(a->mDel, g)]);
      s += a->y[iG(n, g)] * tmp - exp(a->c[iN(a->mC, n)] + 
          a->eps[iNG(a->mEps, n, g)] + tmp);
    }
  }
 
  if(arg * arg > 1e-6){
    tmp = -pow(arg - a->theAlp[a->mTheAlp], 2) / (2 * pow(a->sigAlp[a->mSigAlp], 2)) -
                log(1 - a->piAlp[a->mPiAlp]);
  } else {
    tmp = log(a->piAlp[a->mPiAlp]);
  }

  return s + tmp;
}

__global__ void sampleAlp_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
 
  int g = IDX, G = a->G;
  num_t old, nw, dl, lp, lu;

  if(g < G){ 
    old = a->alp[iG(a->mAlp, g)];
    nw = alpProp(a, g);
    
    dl = lAlp(a, g, nw) - lAlp(a, g, old);
    lp = 0 < dl ? 0 : dl;
    lu = log(runiformDevice(a, g, 0, 1));
    
    if(lu < lp){ /* accept */
      a->alp[iG(a->mAlp + 1, g)] = nw;
      
      if(a->mAlp >= a->burnin)
        ++a->accAlp[g];
    } else { /* reject */
      a->alp[iG(a->mAlp + 1, g)] = old;
    }
  }
} 

__global__ void sampleAlp_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  ++a->mAlp;
}

void sampleAlp(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */

  float myTime;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  fprintf(cfg->log, "alp ");
/*
  sampleAlp_kernel1<<<G_GRID, G_BLOCK>>>(dev_a);
  sampleAlp_kernel2<<<1, 1>>>(dev_a);
  */
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  
  fprintf(cfg->time, "%0.3f ", myTime/MILLISECS); /* elapsed time */  
}