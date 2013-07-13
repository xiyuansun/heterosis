#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void sampleSigAlp_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g, G = a->G;

  for(g = 0; g < a->G; ++g){
    if(pow(a->alp[iG(a->mAlp, g)], 2) > 1e-6){
      a->tmp1[g] = pow(a->alp[iG(a->mAlp, g)] - a->theAlp[a->mTheAlp], 2);
      a->tmp2[g] = 1;
    } else {
      a->tmp1[g] = 0;
      a->tmp2[g] = 0;
    }
  }
}

void sampleSigAlp_kernel2(Chain *a){ /* parallel pairwise sum in Thrust */
  int g;
  num_t rate = 0;  

  for(g = 0; g < a->G; ++g)
    rate += a->tmp1[g];

  a->s1 = rate;
}

void sampleSigAlp_kernel3(Chain *a){ /* parallel pairwise sum in Thrust */
  int g, Galp = 0;
  
  for(g = 0; g < a->G; ++g) 
    Galp += a->tmp2[g];  

  a->s2 = Galp;
}

void sampleSigAlp_kernel4(Chain *a){ /* parallel pairwise sum in Thrust */
  num_t shape = (a->s2 - 1) / 2;
  num_t rate = a->s1 / 2;
  num_t lb = 1/pow(a->sigAlp0, 2);

  if(shape >= 1 && rate > 0){
    a->sigAlp[a->mSigAlp + 1] = 1/sqrt(rgamma(shape, rate, lb));
  } else {
    a->sigAlp[a->mSigAlp + 1] = a->sigAlp[a->mSigAlp]; 
  }

  ++a->mSigAlp;
}

void sampleSigAlp(Chain *a, Config *cfg){ /* host */
  double time;
  clock_t start = clock();

  if(cfg->constSigAlp)
    return;

  if(cfg->verbose)
    printf("sigAlp ");

  sampleSigAlp_kernel1(a);
  sampleSigAlp_kernel2(a);
  sampleSigAlp_kernel3(a);
  sampleSigAlp_kernel4(a); 

  time = ((double) clock() - start) / (SECONDS * CLOCKS_PER_SEC);
  fprintf(cfg->time, "%0.3f ", time);
}