#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

num_t alpProp(Chain *a, int g){ /* device */
      
  num_t gam = a->gamAlp;
  num_t sig = a->sigAlp[a->mSigAlp];

  num_t gprec = 1/(gam * gam);
  num_t sprec = 1/(sig * sig);

  num_t avg = (a->alp[a->mAlp][g] * sprec) / (gprec + sprec);
  num_t s = gam * gam + sig * sig;
  num_t u = uniformHost(0, 1);
  num_t new;
  
  if(u < a->piAlp[a->mPiAlp]){
    new = 0;
  } else {
    new = normalHost(avg, s);
  }

  return new;
}

num_t lAlp(Chain *a, int g, num_t arg){ /* device */
  
  int n;
  num_t s = 0, tmp;
   
  for(n = 0; n < a->N; ++n){
    if(a->grp[n] != 2){
      tmp = mu(a, n, a->phi[a->mPhi][g], arg, a->del[a->mDel][g]);
      s = s + a->y[n][g] * tmp - exp(a->c[a->mC][n] + 
          a->eps[a->mEps][n][g] + tmp);
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

void sampleAlp_kernel1(Chain *a){ /* kernel <<<G, 1>>> */

  int g;
  num_t old, new, dl, lp, lu;

  for(g = 0; g < a->G; ++g){ 

    old = a->alp[a->mAlp][g];
    new = alpProp(a, g);
    
    dl = lAlp(a, g, new) - lAlp(a, g, old);
    lp = 0 < dl ? 0 : dl;
    lu = log(uniformHost(0, 1));
    
    if(lu < lp){ /* accept */
      a->alp[a->mAlp + 1][g] = new;
      a->accAlp[g] = a->accAlp[g] + 1;
    } else { /* reject */
      a->alp[a->mAlp + 1][g] = old;
    }
  }
}

void sampleAlp_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  a->mAlp = a->mAlp + 1;
}


void sampleAlp(Chain *a){ /* host */
  sampleAlp_kernel1(a);
  sampleAlp_kernel2(a);
}