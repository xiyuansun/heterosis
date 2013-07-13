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
    if(pow(a->del[g], 2) > 1e-6){
      a->tmp1[g] = 1;
      a->tmp2[g] = a->del[g];
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
  num_t ss = pow(a->sigDel, 2);
  num_t den = a->s1 * gs + ss;

  num_t m = gs * a->s2 / den;
  num_t s = sqrt(gs * ss / den);

  a->theDel = rnormal(m, s);
}

void sampleTheDel(Chain *a, Config *cfg){ /* host */

  clock_t start = clock();

  if(cfg->constTheDel || !cfg->heterosis)
    return;

  if(cfg->verbose)
    printf("theDel ");

  sampleTheDel_kernel1(a);
  sampleTheDel_kernel2(a);
  sampleTheDel_kernel3(a);
  sampleTheDel_kernel4(a);
  
  cfg->timeTheDel = ((double) clock() - start) / (SECONDS * CLOCKS_PER_SEC);
}