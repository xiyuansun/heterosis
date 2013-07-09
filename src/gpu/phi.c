#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

num_t lPhi(Chain *a, int g, num_t arg){ /* device */
  int n, N = a->N, G = a->G;
  num_t ret, s = 0, tmp = 0; 

  for(n = 0; n < a->N; ++n){
    tmp = mu(a, n, arg, a->alp[a->mAlp][g], a->del[a->mDel][g]);
    s += a->y[iG(n, g)] * tmp - exp(a->c[iN(a->mC, n)] + 
        a->eps[iNG(a->mEps, n, g)] + tmp);
  }
 
  ret = s - pow(arg - a->thePhi[a->mThePhi], 2) / (2 * pow(a->sigPhi[a->mSigPhi], 2));
  return ret;
}

void samplePhi_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g;
  num_t old, new, dl, lp, lu;
  
  for(g = 0; g < a->G; ++g){ 

    old = a->phi[a->mPhi][g];
    new = rnormal(old, a->tunePhi[g]);

    dl = lPhi(a, g, new) - lPhi(a, g, old);
    lp = 0 < dl ? 0 : dl;
    lu = log(runiform(0, 1));
    
    if(lu < lp){ /* accept */
      a->phi[a->mPhi + 1][g] = new;
      a->tunePhi[g] *= 1.1; 
      
      if(a->mPhi >= a->burnin)
        ++a->accPhi[g];
    } else { /* reject */
      a->phi[a->mPhi + 1][g] = old;
      a->tunePhi[g] /= 1.1; 
    }
  }
}

void samplePhi_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  ++a->mPhi;
}

void samplePhi(Chain *a){ /* host */
  samplePhi_kernel1(a);
  samplePhi_kernel2(a);
}
