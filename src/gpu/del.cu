#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

__device__ num_t lDel(Chain *a, int g, num_t arg){ /* device */ 
  int n, N = a->N, G = a->G;
  num_t s = 0, tmp; 
  
  for(n = 0; n < a->N; ++n){
    if(a->grp[n] != 2){
      tmp = mu(a, n, a->phi[iG(a->mPhi, g)], a->alp[iG(a->mAlp, g)], arg);
      s += a->y[iG(n, g)] * tmp - exp(a->c[iN(a->mC, n)] + 
          a->eps[iNG(a->mEps, n, g)] + tmp);
    }
  }
 
  if(arg * arg > 1e-6){
    tmp = -pow(arg - a->theDel[a->mTheDel], 2) / (2 * pow(a->sigDel[a->mSigDel], 2)) -
                log(1 - a->piDel[a->mPiDel]);
  } else {
    tmp = log(a->piDel[a->mPiDel]);
  }

  return s + tmp;
}

__global__ void sampleDel_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g = IDX, G = a->G;
  num_t old, nw, dl, lp, lu;

  if(g < G){ 

    old = a->del[iG(a->mDel, g)];
    nw = delProp(a, g);
    
    dl = lDel(a, g, nw) - lDel(a, g, old);
    lp = 0 < dl? 0 : dl;
    lu = log(runiformDevice(a, g, 0, 1));
    
    if(lu < lp){ /* accept */
      a->del[iG(a->mDel + 1, g)] = nw;
      
      if(a->mDel >= a->burnin)
        ++a->accDel[g];
    } else { /* reject */
      a->del[iG(a->mDel + 1, g)] = old;
    }
  }
}

__global__ void sampleDel_kernel2(Chain *a){ /* kernel <<<1 1>>> */
  ++a->mDel;
}

__host__ void sampleDel(Chain *host_a, Chain *dev_a, Config* cfg){ /* host */
  fprintf(cfg->log, "del ");

  sampleDel_kernel1<<<G_GRID, G_BLOCK>>>(dev_a);
  sampleDel_kernel2<<<1, 1>>>(dev_a);
}