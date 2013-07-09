#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

num_t delProp(Chain *a, int g){ /* device */
      
  num_t gam = a->gamDel;
  num_t sig = a->sigDel[a->mSigDel];

  num_t gprec = 1/(gam * gam);
  num_t sprec = 1/(sig * sig);

  num_t avg = (a->del[a->mDel][g] * sprec) / (gprec + sprec);
  num_t s = gam * gam + sig * sig;
  num_t u = runiform(0, 1);
  num_t new;

  if(u < a->piDel[a->mPiDel]){
    new = 0;
  } else {
    new = rnormal(avg, s);
  }

  return new;
}

num_t lDel(Chain *a, int g, num_t arg){ /* device */ 
  int n, G = a->G;
  num_t s = 0, tmp; 
  
  for(n = 0; n < a->N; ++n){
    if(a->grp[n] != 2){
      tmp = mu(a, n, a->phi[a->mPhi][g], a->alp[a->mAlp][g], arg);
      s += a->y[iNG(n, g)] * tmp - exp(a->c[a->mC][n] + 
          a->eps[a->mEps][n][g] + tmp);
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

void sampleDel_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g;
  num_t old, new, dl, lp, lu;

  for(g = 0; g < a->G; ++g){ 

    old = a->del[a->mDel][g];
    new = delProp(a, g);
    
    dl = lDel(a, g, new) - lDel(a, g, old);
    lp = 0 < dl? 0 : dl;
    lu = log(runiform(0, 1));
    
    if(lu < lp){ /* accept */
      a->del[a->mDel + 1][g] = new;
      
      if(a->mDel >= a->burnin)
        ++a->accDel[g];
    } else { /* reject */
      a->del[a->mDel + 1][g] = old;
    }
  }
}

void sampleDel_kernel2(Chain *a){ /* kernel <<<1 1>>> */
  ++a->mDel;
}

void sampleDel(Chain *a){ /* host */
  sampleDel_kernel1(a);
  sampleDel_kernel2(a);
}