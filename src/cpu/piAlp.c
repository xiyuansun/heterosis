#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void samplePiAlp_kernel1(Chain *a){ /* kernel <<<1, 1>>> */
  int g, G = a->G;

  for(g = 0; g < a->G; ++g){ 
    if(pow(a->alp[iG(a->mAlp, g)], 2) > 1e-6){
      a->tmp1[g] = 1;
    } else {
      a->tmp1[g] = 0;
    }
  }
}

void samplePiAlp_kernel2(Chain *a){ /* pairwise sum in Thrust */
  int g, Galp = 0;
  
  for(g = 0; g < a->G; ++g)
    Galp += a->tmp1[g];

  a->s1 = Galp; 
}

void samplePiAlp_kernel3(Chain *a){ /* kernel <<<1, 1>>> */
  a->piAlp[a->mPiAlp + 1] = rbeta(a->G + a->s1 + a->aTau, a->s1 + a->bTau);
  ++a->mPiAlp;
}

void samplePiAlp(Chain *a, Config *cfg){ /* host */
  fprintf(cfg->log, "piAlp ");

  if(cfg->constPiAlp)
    return;

  samplePiAlp_kernel1(a);
  samplePiAlp_kernel2(a);  
  samplePiAlp_kernel3(a);
}