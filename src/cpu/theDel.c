#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void sampleTheDel_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g, G = a->G;

  for(g = 0; g < a->G; ++g){ 
    if(pow(a->del[iG(a->mDel, g)], 2) > 1e-6){
      a->tmp1[g] = 1;
      a->tmp2[g] = a->del[iG(a->mDel, g)];
    } else {
      a->tmp1[g] = 0;
      a->tmp2[g] = 0;
    }
  }
}

void sampleTheDel_kernel2(Chain *a){ /* pairwise sum in Thrust */
  int g, Gdel = 0;
  
  for(g = 0; g < a->G; ++g)   
    Gdel += a->tmp1[g];

  a->s1 = Gdel;
}

void sampleTheDel_kernel3(Chain *a){ /* pairwise sum in Thrust */
  int g;
  num_t sm = 0;
  
  for(g = 0; g < a->G; ++g) 
    sm += a->tmp2[g];

  a->s2 = sm;
}

void sampleTheDel_kernel4(Chain *a){ /* kernel <<<1, 1>>> */

  num_t gs = pow(a->gamDel, 2);
  num_t ss = pow(a->sigDel[a->mSigDel], 2);
  num_t den = a->s1 * gs + ss;

  num_t m = gs * a->s2 / den;
  num_t s = sqrt(gs * ss / den);

  a->theDel[a->mTheDel + 1] = rnormal(m, s);
  ++a->mTheDel;
}

void sampleTheDel(Chain *a, Config *cfg){ /* host */
  double time;
  clock_t start = clock();

  if(cfg->constTheDel || !cfg->heterosis)
    return;

  fprintf(cfg->log, "theDel ");

  sampleTheDel_kernel1(a);
  sampleTheDel_kernel2(a);
  sampleTheDel_kernel3(a);
  sampleTheDel_kernel4(a);
  
  time = ((double) clock() - start) / (SECS * CLOCKS_PER_SEC);
  fprintf(cfg->time, "%0.3f ", time);  
}