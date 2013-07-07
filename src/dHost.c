#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void lD_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g;

  for(g = 0; g < a->G; ++g){ 
    a->tmp1[g] = 2 * log(a->eta[a->mEta][g]);
    a->tmp2[g] = 1/(a->eta[a->mEta][g] * a->eta[a->mEta][g]);
  }
}

void lD_kernel2(Chain *a){ /* kernel: pairwise sum in Thrust */
  int g;
  a->s1 = 0;

  for(g = 0; g < a->G; ++g) /* PARALLELIZE */
    a->s1 = a->s1 + a->tmp1[g];
}

void lD_kernel3(Chain *a){ /* kernel: pairwise sum in Thrust */
  int g;
  a->s2 = 0;

  for(g = 0; g < a->G; ++g) /* PARALLELIZE */
    a->s2 = a->s2 + a->tmp2[g];
}

void lD_kernel4(Chain *a, int newArg){ /* kernel <<<1, 1>>> */
  num_t arg, ret;
 
  if(newArg){
    arg = a->New[0][0];
  } else{
    arg = a->Old[0][0];
  }

  a->tmp1[1] = arg * a->tau[a->mTau] * a->tau[a->mTau] / 2;

  ret = -a->G * lgamma(arg/2) + (a->G * arg / 2) * log(a->tmp1[1]);
  ret = ret  - (arg/2 + 1) * a->s1 - a->tmp1[1] * a->s2;

  if(newArg){
    a->lNew[0][0] = ret;
  } else{
    a->lOld[0][0] = ret;
  }
}

void lD(Chain *a, int newArg){ /* host */
  
  if(newArg){
    if(a->New[0][0] <= 0 || a->New[0][0] > a->d0)
      a->lNew[0][0] = NUM_TMIN;
  } else {
    if(a->Old[0][0] <= 0 || a->Old[0][0] > a->d0)
      a->lOld[0][0] = NUM_TMIN; 
  }

  lD_kernel1(a);
  lD_kernel2(a);
  lD_kernel3(a);
  lD_kernel4(a, newArg);
}

void sampleD_kernel1(Chain *a){ /* kernel <<<1, 1>>> */
  a->Old[0][0] = a->d[a->mD];
  
  do {
    a->New[0][0] = normalHost(a->Old[0][0], a->tuneD);
  } while(a->New[0][0] < 0);
}

void sampleD_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  num_t dl = a->lNew[0][0] - a->lOld[0][0];
  num_t lp = 0 < dl ? 0 : dl;
  num_t lu = log(uniformHost(0, 1));

  if(lu < lp){ /* accept */
    a->d[a->mD + 1] = a->New[0][0];
    a->tuneD = a->tuneD * 1.1; /* Increase the proposal variance to avoid getting 
                                  stuck in a mode */
    a->accD = a->accD + 1;
  } else { /* reject */
    a->d[a->mD + 1] = a->Old[0][0];
    a->tuneD = a->tuneD / 1.1; /* If you're rejecting too often, decrease the proposal 
                                  variance to sample closer to the last accepted value. */
  }

  a->mD = a->mD + 1;
}

void sampleD(Chain *a, Config *cfg){ /* host */
  if(cfg->constD)
    return;
   
  sampleD_kernel1(a);

  lD(a, 1);
  lD(a, 0);

  sampleD_kernel2(a);
}
