#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

num_t lEps(Chain *a, int n, int g, num_t arg){ /* device */
  int M = a->M, N = a->N, G = a->G;

  return a->y[iNG(n, g)] * arg - exp(a->c[iMN(a->mC, n)] + arg + mu(a, n, a->phi[iMG(a->mPhi, g)], 
                                     a->alp[iMG(a->mAlp, g)], a->del[iMG(a->mDel, g)])) 
                          - (arg * arg) / (2 * a->eta[iMG(a->mEta, g)] * a->eta[iMG(a->mEta, g)]);
}

void sampleEps_kernel1(Chain *a){ /* kernel <<<N, G>>> */
  int n, g;
  int M = a->M, N = a->N, G = a->G;
  num_t Old;
  num_t New;
  num_t dl;
  num_t lp;
  num_t lu;

  for(g = 0; g < G; ++g){
    for(n = 0; n < N; ++n){ 
      Old = a->eps[iMNG(a->mEps, n, g)];
      New = rnormal(Old, a->tuneEps[iNG(n, g)]);

      dl = lEps(a, n, g, New) - lEps(a, n, g, Old);
      lp = 0 < dl ? 0 : dl;
      lu = log(runiform(0, 1));
      
      if(lu < lp){ /* accept */
        a->eps[iMNG(a->mEps + 1, n, g)] = New;
        a->tuneEps[iNG(n, g)] = a->tuneEps[iNG(n, g)] * 1.1;
        
        if(a->mEps >= a->burnin)
          ++a->accEps[iNG(n, g)]; 
      } else { /* reject */
        a->eps[iMNG(a->mEps + 1, n, g)] = Old;
        a->tuneEps[iNG(n, g)] = a->tuneEps[iNG(n, g)] / 1.1;
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