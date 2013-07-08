#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void sampleSigC(Chain *a){ /* kernel <<<1, 1>>> */
  int n;
  int M = a->M, N = a->N, G = a->G;
  num_t rate, shape, lb;

  if(a->constSigC)
    return; 

  rate = 0;
  for(n = 0; n < a->N; ++n) 
    rate += a->c[a->mC][n] * a->c[a->mC][n];
  
  shape = (a->N - 1) / 2; 
  rate = rate / 2;
  lb = 1 / pow(a->sigC0, 2); 

  if(shape >= 1 && rate > 0){
    a->sigC[a->mSigC + 1] = 1/sqrt(rgamma(shape, rate, lb));
  } else {
    a->sigC[a->mSigC + 1] = a->sigC[a->mSigC];
  }

  ++a->mSigC;
}