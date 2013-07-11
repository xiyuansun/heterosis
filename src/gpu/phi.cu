#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

__device__ num_t lPhi(Chain *a, int g, num_t arg){ /* device */
  int n, N = a->N, G = a->G;
  num_t ret, s = 0, tmp = 0; 

  for(n = 0; n < a->N; ++n){
    tmp = mu(a, n, arg, a->alp[iG(a->mAlp, g)], a->del[iG(a->mDel, g)]);
    s += a->y[iG(n, g)] * tmp - exp(a->c[iN(a->mC, n)] + 
        a->eps[iNG(a->mEps, n, g)] + tmp);
  }
 
  ret = s - pow(arg - a->thePhi[a->mThePhi], 2) / (2 * pow(a->sigPhi[a->mSigPhi], 2));
  return ret;
}

__global__ void samplePhi_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g = GENE, G = a->G;
  num_t old, nw, dl, lp, lu;
  
  if(g < G){ 
    old = a->phi[iG(a->mPhi, g)];
    nw = rnormalDevice(a, g, old, a->tunePhi[g]);

    dl = lPhi(a, g, nw) - lPhi(a, g, old);
    lp = 0 < dl ? 0 : dl;
    lu = log(runiformDevice(a, g, 0, 1));
    
    if(lu < lp){ /* accept */
      a->phi[iG(a->mPhi + 1, g)] = nw;
      a->tunePhi[g] *= 1.1; 
      
      if(a->mPhi >= a->burnin)
        ++a->accPhi[g];
    } else { /* reject */
      a->phi[iG(a->mPhi + 1, g)] = old;
      a->tunePhi[g] /= 1.1; 
    }
  }
}

__global__ void samplePhi_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  ++a->mPhi;
}

void samplePhi(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */
  samplePhi_kernel1<<<NBLOCKS, NTHREADS>>>(dev_a);
  samplePhi_kernel2<<<1, 1>>>(dev_a);
}
