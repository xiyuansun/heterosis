#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void sampleThePhi_kernel1(Chain *a){ /* pairwise sum in Thrust */
  int g;
  
  a->s1 = 0; 
  for(g = 0; g < a->G; ++g)
    a->s1 = a->s1 + a->phi[a->mPhi][g];
}

void sampleThePhi_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  num_t gs = a->gamPhi * a->gamPhi;
  num_t ss = a->sigPhi[a->mSigPhi] * a->sigPhi[a->mSigPhi];
  num_t den = (a->G * gs + ss);

  num_t m = gs * a->s1 / den;
  num_t s = gs * ss / den;

  a->thePhi[a->mThePhi + 1] = normalHost(m, s);
  a->mThePhi = a->mThePhi + 1;
}

void sampleThePhi(Chain *a){ /* host */
  sampleThePhi_kernel1(a);
  sampleThePhi_kernel2(a);
}