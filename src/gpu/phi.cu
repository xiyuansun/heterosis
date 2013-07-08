#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

num_t lPhi(Chain *a, int g, num_t arg){ /* device */
  int n;
  int N = a->N, G = a->G;
  
  num_t ret, s = 0, tmp = 0; 

  for(n = 0; n < N; ++n){
    tmp = mu(a, n, arg, a->alp[iMG(a->mAlp, g)], a->del[iMG(a->mDel, g)]);
    s += a->y[iNG(n, g)] * tmp - exp(a->c[iMG(a->mC, n)] + 
        a->eps[iMNG(a->mEps, n, g)] + tmp);
  }
 
  ret = s - pow(arg - a->thePhi[a->mThePhi], 2) / (2 * pow(a->sigPhi[a->mSigPhi], 2));
  return ret;
}

void samplePhi_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g;
  int G = a->G;
  
  num_t Old, New, dl, lp, lu;
  
  for(g = 0; g < G; ++g){ 

    Old = a->phi[iMG(a->mPhi, g)];
    New = rnormal(Old, a->tunePhi[g]);

    dl = lPhi(a, g, New) - lPhi(a, g, Old);
    lp = 0 < dl ? 0 : dl;
    lu = log(runiform(0, 1));
    
    if(lu < lp){ /* accept */
      a->phi[iMG(a->mPhi + 1, g)] = New;
      a->tunePhi[g] *= 1.1; 
      
      if(a->mPhi >= a->burnin)
        ++a->accPhi[g];
    } else { /* reject */
      a->phi[iMG(a->mPhi + 1, g)] = Old;
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
