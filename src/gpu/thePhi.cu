#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void sampleThePhi_kernel1(Chain *a){ /* pairwise sum in Thrust */
  int g, G = a->G;
  
  a->tmp1[0] = 0; 
  for(g = 0; g < a->G; ++g)
    a->tmp1[0] += a->phi[iG(a->mPhi, g)];
}

void sampleThePhi_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  num_t gs = a->gamPhi * a->gamPhi;
  num_t ss = a->sigPhi[a->mSigPhi] * a->sigPhi[a->mSigPhi];
  num_t den = (a->G * gs + ss);

  num_t m = gs * a->tmp1[0] / den;
  num_t s = gs * ss / den;

  a->thePhi[a->mThePhi + 1] = rnormal(m, s);
  ++a->mThePhi;
}

void sampleThePhi(Chain *a, Config *cfg){ /* host */
  if(cfg->constThePhi)
    return;

  sampleThePhi_kernel1(a);
  sampleThePhi_kernel2(a);
}