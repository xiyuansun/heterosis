#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

__device__ num_t lEps(Chain *a, int n, int g, num_t arg){ /* device */
  int N = a->N, G = a->G;
  return a->y[iG(n, g)] * arg - exp(a->c[iN(a->mC, n)] + arg + mu(a, n, a->phi[iG(a->mPhi, g)], 
                                     a->alp[iG(a->mAlp, g)], a->del[iG(a->mDel, g)])) 
                          - (arg * arg) / (2 * pow(a->eta[iG(a->mEta, g)], 2));
}

__global__ void sampleEps_kernel1(Chain *a){ /* kernel <<<N, G>>> */
  int n = IDY, g = IDX, N = a->N, G = a->G, id = ID;
  num_t old, nw, dl, lp, lu;

  if(g < G){
    if(n < N){ 
      old = a->eps[iNG(a->mEps, n, g)];
      nw = rnormalDevice(a, id, old, a->tuneEps[iG(n, g)]);

      dl = lEps(a, n, g, nw) - lEps(a, n, g, old);
      lp = 0 < dl ? 0 : dl;
      lu = log(runiformDevice(a, id, 0, 1));
      
      if(lu < lp){ /* accept */
        a->eps[iNG(a->mEps + 1, n, g)] = nw;
        a->tuneEps[iG(n, g)] *= 1.1;
        
        if(a->mEps >= a->burnin)
          ++a->accEps[iG(n, g)]; 
      } else { /* reject */
        a->eps[iNG(a->mEps + 1, n, g)] = old;
        a->tuneEps[iG(n, g)] /= 1.1;
      }
    }
  }
}

__global__ void sampleEps_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  ++a->mEps;
}

void sampleEps(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */
  sampleEps_kernel1<<<GN_GRID, GN_BLOCK>>>(dev_a);
  sampleEps_kernel2<<<1, 1>>>(dev_a);
}