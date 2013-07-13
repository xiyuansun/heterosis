#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void sampleEta_kernel1(Chain *a){ /* kernel <<<1, 1>>> */
  a->s1 = (a->N + a->d[a->mD]) / 2; 
}

void sampleEta_kernel2(Chain *a){ /* kernel <<<G, 1>>> */
  int n, g, N = a->N, G = a->G;
  num_t shape = a->s1, rate;

  for(g = 0; g < a->G; ++g){

    rate = 0;
    for(n = 0; n < a->N; ++n) 
      rate += pow(a->eps[iNG(a->mEps, n, g)], 2);
  
    rate = (rate + a->d[a->mD] * a->tau[a->mTau] * a->tau[a->mTau]) / 2; 

    if(shape >= 1 && rate > 0){
      a->eta[iG(a->mEta + 1, g)] = 1/sqrt(rgamma(shape, rate, 0));
    } else {
      a->eta[iG(a->mEta + 1, g)] = a->eta[iG(a->mEta, g)];
    }
  }
}

void sampleEta_kernel3(Chain *a){ /* kernel <<<1, 1>>> */
  ++a->mEta;
}

void sampleEta(Chain *a, Config *cfg){
  double time;
  clock_t start = clock();

  if(cfg->verbose)
    printf("eta ");

  sampleEta_kernel1(a);
  sampleEta_kernel2(a);
  sampleEta_kernel3(a);

  time = ((double) clock() - start) / (SECONDS * CLOCKS_PER_SEC);
  fprintf(cfg->time, "%0.3f ", time);
}