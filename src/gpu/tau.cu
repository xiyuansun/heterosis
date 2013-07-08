#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void sampleTau_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g;
  int M = a->M, N = a->N, G = a->G;
  
  for(g = 0; g < a->G; ++g)
    a->tmp1[g] = 1/pow(a->eta[iMG(a->mEta, g)], 2);
}

void sampleTau_kernel2(Chain *a){ /* pairwise sum in Thrust */
  int g;
  int M = a->M, N = a->N, G = a->G;
  num_t tmp = 0;
  
  for(g = 0; g < a->G; ++g) 
    tmp += a->tmp1[g];

  a->s1 = tmp;
}

void sampleTau_kernel3(Chain *a){ /* kernel<<<1, 1>>> */
  num_t rate = a->s1 * a->d[a->mD] / 2 + a->bTau;
  num_t shape = a->aTau + a->G * a->d[a->mD] / 2;

  if(shape >= 1 && rate > 0){
    a->tau[a->mTau + 1] = 1/sqrt(rgamma(shape, rate, 0));
  } else {
    a->tau[a->mTau + 1] = a->tau[a->mTau];
  }

  ++a->mTau;
}

void sampleTau(Chain *a, Config *cfg){ /* host */
  if(cfg->constTau)
    return;

  sampleTau_kernel1(a);
  sampleTau_kernel2(a);
  sampleTau_kernel3(a);
}