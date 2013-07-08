#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void sampleSigPhi_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g;
  int G = a->G;

  for(g = 0; g < a->G; ++g) 
    a->tmp1[g] = pow(a->phi[iMG(a->mPhi, g)] - a->thePhi[a->mThePhi], 2);
}

void sampleSigPhi_kernel2(Chain *a){ /* parallel pairwise sum in Thrust */
  int g;
  int G = a->G;
  
  num_t rate = 0;
  
  for(g = 0; g < a->G; ++g)
    rate += a->tmp1[g];
  a->s1 = rate;  
}

void sampleSigPhi_kernel3(Chain *a){ /* kernel <<<1, 1>>> */
  num_t rate = a->s1 / 2;
  num_t shape = (a->G - 1) / 2;
  num_t lb = 1/pow(a->sigPhi0, 2);

  if(shape >= 1 && rate > 0){
    a->sigPhi[a->mSigPhi + 1] = 1/sqrt(rgamma(shape, rate, lb));
  } else {
    a->sigPhi[a->mSigPhi + 1] = a->sigPhi[a->mSigPhi];
  }

  ++a->mSigPhi;
}

void sampleSigPhi(Chain *a, Config *cfg){ /* host */
  if(cfg->constSigPhi)
    return;

  sampleSigPhi_kernel1(a);
  sampleSigPhi_kernel2(a);
  sampleSigPhi_kernel3(a);
}