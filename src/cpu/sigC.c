#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void sampleSigC(Chain *a, Config *cfg){ /* kernel <<<1, 1>>> */

  int n;
  num_t rate, shape, lb;
  clock_t start = clock();

  if(a->constSigC)
    return; 

  if(cfg->verbose)
    printf("sigC ");

  rate = 0;
  for(n = 0; n < a->N; ++n) 
    rate += a->c[n] * a->c[n];
  
  shape = (a->N - 1) / 2; 
  rate = rate / 2;
  lb = 1 / pow(a->sigC0, 2); 

  if(shape > 0 && rate > 0){
    a->sigC = 1/sqrt(rgamma(shape, rate, lb));
  } 

  cfg->timeSigC = ((num_t) clock() - start) / (SECONDS * CLOCKS_PER_SEC);
}
