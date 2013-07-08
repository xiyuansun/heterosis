#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

num_t lEps(Chain *a, int n, int g, num_t arg){ /* device */

  return a->y[n][g] * arg - exp(a->c[a->mC][n] + arg + mu(a, n, a->phi[a->mPhi][g], 
                                     a->alp[a->mAlp][g], a->del[a->mDel][g])) 
                          - (arg * arg) / (2 * a->eta[a->mEta][g] * a->eta[a->mEta][g]);
}

void sampleEps_kernel1(Chain *a){ /* kernel <<<N, G>>> */
  int n, g;
  num_t old, new, dl, lp, lu;

  for(g = 0; g < a->G; ++g){
    for(n = 0; n < a->N; ++n){ 
      old = a->eps[a->mEps][n][g];
      new = normalHost(old, a->tuneEps[n][g]);

      dl = lEps(a, n, g, new) - lEps(a, n, g, old);
      lp = 0 < dl ? 0 : dl;
      lu = log(uniformHost(0, 1));
      
      if(lu < lp){ /* accept */
        a->eps[a->mEps + 1][n][g] = new;
        a->tuneEps[n][g] = a->tuneEps[n][g] * 1.1;
        
        if(a->mEps >= a->burnin)
          ++a->accEps[n][g]; 
      } else { /* reject */
        a->eps[a->mEps + 1][n][g] = old;
        a->tuneEps[n][g] = a->tuneEps[n][g] / 1.1;
      }
    }
  }
}

void sampleEps_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  ++a->mEps;
}

void sampleEps(Chain *a){ /* host */
  sampleEps_kernel1(a);
  sampleEps_kernel2(a);
}