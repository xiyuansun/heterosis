#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void sampleThePhi_kernel1(Chain *a){ /* pairwise sum in Thrust */
  int g;
  
  a->tmp1[0] = 0; 
  for(g = 0; g < a->G; ++g)
    a->tmp1[0] += a->phi[g];
}

void sampleThePhi_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  num_t gs = a->gamPhi * a->gamPhi;
  num_t ss = a->sigPhi * a->sigPhi;
  num_t den = (a->G * gs + ss);

  num_t m = gs * a->tmp1[0] / den;
  num_t s = sqrt(gs * ss / den);

  a->thePhi = rnormal(m, s);
}

void sampleThePhi(Chain *a, Config *cfg){ /* host */  

  clock_t start = clock();

  if(cfg->constThePhi)
    return;

  if(cfg->verbose)
    printf("thePhi ");

  if(!cfg->phiPrior){
	sampleThePhi_kernel1(a);
	sampleThePhi_kernel2(a);
  }
  
  cfg->timeThePhi = ((double) clock() - start) / (SECONDS * CLOCKS_PER_SEC);
}
