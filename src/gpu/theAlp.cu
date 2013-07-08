#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void sampleTheAlp_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g;
  int M = a->M, N = a->N, G = a->G;

  for(g = 0; g < a->G; ++g){
    if(pow(a->alp[a->mAlp][g], 2) > 1e-6){
      a->tmp1[g] = 1;
      a->tmp2[g] = a->alp[a->mAlp][g];
    } else {
      a->tmp1[g] = 0;
      a->tmp2[g] = 0;
    }
  }
}

void sampleTheAlp_kernel2(Chain *a){ /* parallel pairwise sum in Thrust */
  int g, Galp = 0;
  int M = a->M, N = a->N, G = a->G;
  
  for(g = 0; g < a->G; ++g) 
    Galp += a->tmp1[g];

  a->s1 = Galp;
}

void sampleTheAlp_kernel3(Chain *a){ /* parallel pairwise sum in Thrust */
  int g;
  int M = a->M, N = a->N, G = a->G;
  num_t sm = 0;
  
  for(g = 0; g < a->G; ++g) 
    sm += a->tmp2[g];

  a->s2 = sm;
}

void sampleTheAlp_kernel4(Chain *a){ /* kernel <<<1, 1>>> */

  num_t gs = pow(a->gamAlp, 2);
  num_t ss = pow(a->sigAlp[a->mSigAlp], 2);
  num_t den = a->s1 * gs + ss;

  num_t m = gs * a->s2 / den;
  num_t s = sqrt(gs * ss / den);

  a->theAlp[a->mTheAlp + 1] = rnormal(m, s);
  ++a->mTheAlp;
}

void sampleTheAlp(Chain *a, Config *cfg){ /* host */
  if(cfg->constTheAlp)
    return;

  sampleTheAlp_kernel1(a);
  sampleTheAlp_kernel2(a);
  sampleTheAlp_kernel3(a);
  sampleTheAlp_kernel4(a);
}
