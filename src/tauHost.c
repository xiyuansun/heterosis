#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void sampleTau_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g;
  
  for(g = 0; g < a->G; ++g)
    a->tmp1[g] = 1/pow(a->eta[a->mEta][g], 2);
}

void sampleTau_kernel2(Chain *a){ /* pairwise sum in Thrust */
  int g;
  num_t tmp = 0;
  
  for(g = 0; g < a->G; ++g) 
    tmp = tmp + a->tmp1[g];

  a->rate = tmp;
}

void sampleTau_kernel3(Chain *a){ /* kernel<<<1, 1>>> */
  num_t rate = a->rate * a->d[a->mD] / 2 + a->bTau;
  num_t shape = a->aTau + a->G * a->d[a->mD] / 2;

  if(shape >= 1 && rate > 0){
    a->tau[a->mTau + 1] = 1/sqrt(gammaHost(shape, rate, 0));
  } else {
    a->tau[a->mTau + 1] = a->tau[a->mTau];
  }

  a->mTau = a->mTau + 1;
}

void sampleTau(Chain *a){ /* host */
  sampleTau_kernel1(a);
  sampleTau_kernel2(a);
  sampleTau_kernel3(a);
}