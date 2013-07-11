#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

__device__ num_t alpProp(Chain *a, int g){ /* device */
  int G = a->G;

  num_t gam = a->gamAlp;
  num_t sig = a->sigAlp[a->mSigAlp];

  num_t gprec = 1/(gam * gam);
  num_t sprec = 1/(sig * sig);

  num_t avg = (a->alp[iG(a->mAlp, g)] * sprec) / (gprec + sprec);
  num_t s = gam * gam + sig * sig;
  num_t u = runiform(0, 1);
  num_t nw;
  
  if(u < a->piAlp[a->mPiAlp]){
    nw = 0;
  } else {
    nw = rnormal(avg, s);
  }

  return nw;
}

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

  int g = GENE, G = a->G;
  num_t old, nw, dl, lp, lu;

  if(g < G){ 

    old = a->alp[iG(a->mAlp, g)];
    nw = alpProp(a, g);
    
    dl = lAlp(a, g, nw) - lAlp(a, g, old);
    lp = 0 < dl ? 0 : dl;
    lu = log(runiform(0, 1));
    
    if(lu < lp){ /* accept */
      a->alp[iG(a->mAlp + 1, g)] = nw;
      
      if(a->mAlp >= a->burnin)
        ++a->accAlp[g];
    } else { /* reject */
      a->alp[iG(a->mAlp + 1, g)] = old;
    }
  }
}

void sampleAlp_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  ++a->mAlp;
}


void sampleAlp(Chain *a){ /* host */
  sampleAlp_kernel1<<<NBLOCKS, NTHREADS>>>(a);
  sampleAlp_kernel2<<<1, 1>>>(a);
}