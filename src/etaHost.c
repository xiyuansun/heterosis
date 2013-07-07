#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void sampleEta_kernel1(Chain *a){ /* kernel <<<1, 1>>> */
  a->shape = (a->N + a->d[a->mD]) / 2; 
}

void sampleEta_kernel2(Chain *a){ /* kernel <<<G, 1>>> */
  int n, g;
  num_t rate;

  for(g = 0; g < a->G; ++g){

    rate = 0;
    for(n = 0; n < a->N; ++n) 
      rate = rate + a->eps[a->mEps][n][g] * a->eps[a->mEps][n][g];
  
    rate = (rate + a->d[a->mD] * a->tau[a->mTau] * a->tau[a->mTau]) / 2; 

    if(a->shape >= 1 && rate > 0){
      a->eta[a->mEta + 1][g] = 1/sqrt(gammaHost(a->shape, rate, 0));
    } else {
      a->eta[a->mEta + 1][g] = a->eta[a->mEta][g];
    }
  }
}

void sampleEta_kernel3(Chain *a){ /* kernel <<<1, 1>>> */
  a->mEta = a->mEta + 1;
}

void sampleEta(Chain *a){
  sampleEta_kernel1(a);
  sampleEta_kernel2(a);
  sampleEta_kernel3(a);
}