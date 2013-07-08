#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void sampleSigDel_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g;

  for(g = 0; g < a->G; ++g){ 
    if(pow(a->del[a->mDel][g], 2) > 1e-6){
      a->tmp1[g] = pow(a->del[a->mDel][g] - a->theDel[a->mTheDel], 2);
      a->tmp2[g] = 1;
    } else {
      a->tmp1[g] = 0;
      a->tmp2[g] = 0;
    }
  }
}

void sampleSigDel_kernel2(Chain *a){ /* pairwise sum in Thrust */
  int g;
  num_t rate = 0;
  
  for(g = 0; g < a->G; ++g) 
    rate += a->tmp1[g];

  a->s1 = rate;
}

void sampleSigDel_kernel3(Chain *a){ /* pairwise sum in Thrust */
  int g, Gdel = 0;
  
  for(g = 0; g < a->G; ++g) 
    Gdel += a->tmp2[g];

  a->s2 = Gdel;
}

void sampleSigDel_kernel4(Chain *a){ /* kernel <<<1, 1>>> */
  num_t shape = (a->s2 - 1) / 2;
  num_t rate = a->s1 / 2;
  num_t lb = 1/pow(a->sigDel0, 2);

  if(shape >= 1 && rate > 0){
    a->sigDel[a->mSigDel + 1] = 1/sqrt(rgamma(shape, rate, lb));
  } else {
    a->sigDel[a->mSigDel + 1] = a->sigDel[a->mSigDel];
  }

  ++a->mSigDel;
}

void sampleSigDel(Chain *a, Config *cfg){ /* host */
  if(cfg->constSigDel || !cfg->heterosis)
    return;

  sampleSigDel_kernel1(a);
  sampleSigDel_kernel2(a);
  sampleSigDel_kernel3(a);
  sampleSigDel_kernel4(a);
}