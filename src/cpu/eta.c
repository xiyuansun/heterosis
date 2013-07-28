#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void sampleEta_kernel1(Chain *a){ /* kernel <<<1, 1>>> */
  a->s1 = (a->N + a->d) / 2; 
}

void sampleEta_kernel2(Chain *a){ /* kernel <<<G, 1>>> */
  int n, g, G = a->G;
  num_t shape = a->s1, rate;

  for(g = 0; g < a->G; ++g){

    rate = 0;
    for(n = 0; n < a->N; ++n) 
      rate += pow(a->eps[iG(n, g)], 2);
  
    rate = (rate + a->d * a->tau * a->tau) / 2; 

    if(shape > 0 && rate > 0){
      a->eta[g] = 1/sqrt(rgamma(shape, rate, 0));
    } else {
      a->eta[g] = a->eta[g];
    }
  }
}

void sampleEta(Chain *a, Config *cfg){

  clock_t start = clock();

  if(cfg->verbose)
    printf("eta ");

  sampleEta_kernel1(a);
  sampleEta_kernel2(a);

  cfg->timeEta = ((num_t) clock() - start) / (SECONDS * CLOCKS_PER_SEC);
}