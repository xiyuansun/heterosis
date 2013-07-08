#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void sampleEta_kernel1(Chain *a){ /* kernel <<<1, 1>>> */
  a->s1 = (a->N + a->d[a->mD]) / 2; 
}

void sampleEta_kernel2(Chain *a){ /* kernel <<<G, 1>>> */
  int n, g;
  int M = a->M, N = a->N, G = a->G;
  
  num_t shape = a->s1, rate;

  for(g = 0; g < G; ++g){

    rate = 0;
    for(n = 0; n < N; ++n) 
      rate += a->eps[a->mEps][n][g] * a->eps[a->mEps][n][g];
  
    rate = (rate + a->d[a->mD] * a->tau[a->mTau] * a->tau[a->mTau]) / 2; 

    if(shape >= 1 && rate > 0){
      a->eta[a->mEta + 1][g] = 1/sqrt(rgamma(shape, rate, 0));
    } else {
      a->eta[a->mEta + 1][g] = a->eta[a->mEta][g];
    }
  }
}

void sampleEta_kernel3(Chain *a){ /* kernel <<<1, 1>>> */
  ++a->mEta;
}

void sampleEta(Chain *a){
  sampleEta_kernel1(a);
  sampleEta_kernel2(a);
  sampleEta_kernel3(a);
}