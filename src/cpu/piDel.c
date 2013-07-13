#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void samplePiDel_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g, G = a->G;

  for(g = 0; g < a->G; ++g){ 
    if(pow(a->del[iG(a->mDel, g)], 2) > 1e-6){
      a->tmp1[g] = 1; 
    } else {
      a->tmp1[g] = 0;
    }
  } 
}

void samplePiDel_kernel2(Chain *a){ /* pairwise sum in Thrust */

  int g, Gdel = 0;
  for(g = 0; g < a->G; ++g) 
     Gdel += a->tmp1[g];

  a->s1 = Gdel;
}

void samplePiDel_kernel3(Chain *a){ /* kernel <<<1, 1>>> */
  a->piDel[a->mPiDel + 1] = rbeta(a->G + a->s1 + a->aTau, a->s1 + a->bTau);
  ++a->mPiDel;
}

void samplePiDel(Chain *a, Config *cfg){ /* host */
  double time;
  clock_t start = clock();

  if(cfg->constPiDel || !cfg->heterosis)
    return;

  if(cfg->verbose)
    printf("piDel ");

  samplePiDel_kernel1(a);
  samplePiDel_kernel2(a);
  samplePiDel_kernel3(a);

  time = ((double) clock() - start) / (SECS * CLOCKS_PER_SEC);
  fprintf(cfg->time, "%0.3f ", time);
}
