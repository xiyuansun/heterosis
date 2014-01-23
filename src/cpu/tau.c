#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void sampleTau_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g;
  
  for(g = 0; g < a->G; ++g)
    a->tmp1[g] = 1/pow(a->eta[g], 2);
}

void sampleTau_kernel2(Chain *a){ /* pairwise sum in Thrust */
  int g;
  num_t tmp = 0;
  
  for(g = 0; g < a->G; ++g) 
    tmp += a->tmp1[g];

  a->s1 = tmp;
}

void sampleTau_kernel3(Chain *a){ /* kernel<<<1, 1>>> */
  num_t rate = a->s1 * a->d / 2 + a->bTau;
  num_t shape = a->aTau + a->G * a->d / 2;

  if(shape > 0 && rate > 0){
    a->tau = 1/sqrt(rgamma(shape, rate, 0));
  } 
}

void sampleTau(Chain *a, Config *cfg){ /* host */

  clock_t start = clock();

  if(cfg->constTau)
    return;

  if(cfg->verbose)
    printf("tau ");

  sampleTau_kernel1(a);
  sampleTau_kernel2(a);
  sampleTau_kernel3(a);

  cfg->timeTau = ((double) clock() - start) / (SECONDS * CLOCKS_PER_SEC);
}
