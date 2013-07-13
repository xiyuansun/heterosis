#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void lD_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g;

  for(g = 0; g < a->G; ++g){ 
    a->tmp1[g] = 2 * log(a->eta[g]);
    a->tmp2[g] = 1/(a->eta[g] * a->eta[g]);
  }
}

void lD_kernel2(Chain *a){ /* kernel: pairwise sum in Thrust */
  int g;
  a->s1 = 0;

  for(g = 0; g < a->G; ++g) /* PARALLELIZE */
    a->s1 += a->tmp1[g];
}

void lD_kernel3(Chain *a){ /* kernel: pairwise sum in Thrust */
  int g;
  a->s2 = 0;

  for(g = 0; g < a->G; ++g) /* PARALLELIZE */
    a->s2 += a->tmp2[g];
}

void lD_kernel4(Chain *a, int newArg){ /* kernel <<<1, 1>>> */
  num_t arg, ret, tmp;
 
  if(newArg){
    arg = a->New[0];
  } else{
    arg = a->Old[0];
  }

  tmp = arg * a->tau * a->tau / 2;
  ret = -a->G * lgamma(arg/2) + (a->G * arg / 2) * log(tmp);
  ret -= (arg/2 + 1) * a->s1 - tmp * a->s2;
  
  if(ret < 1e-6 || ret > a->d0)
    ret = NUM_TMIN;

  if(newArg){
    a->lNew[0] = ret;
  } else{
    a->lOld[0] = ret;
  }
}

void lD(Chain *a, int newArg){ /* host */

  lD_kernel1(a);
  lD_kernel2(a);
  lD_kernel3(a);
  lD_kernel4(a, newArg);
}

void sampleD_kernel1(Chain *a){ /* kernel <<<1, 1>>> */
  a->Old[0] = a->d;
  
  do {
    a->New[0] = rnormal(a->Old[0], a->tuneD);
  } while(a->New[0] < 1e-6);
}

void sampleD_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  num_t dl = a->lNew[0] - a->lOld[0];
  num_t lp = 0 < dl ? 0 : dl;
  num_t lu = log(runiform(0, 1));

  if(lu < lp){ /* accept */
    a->d = a->New[0];
    a->tuneD *= 1.1; /* Increase the proposal variance to avoid getting 
                                  stuck in a mode */
    
    if(a->m > a->burnin) 
      ++a->accD;
  } else { /* reject */
    a->tuneD /= 1.1; /* If you're rejecting too often, decrease the proposal 
                                  variance to sample closer to the last accepted value. */
  }
}

void sampleD(Chain *a, Config *cfg){ /* host */

  clock_t start = clock();
  
  if(cfg->constD)
    return;
  
  if(cfg->verbose)  
    printf("d ");
   
  sampleD_kernel1(a);

  lD(a, 1);
  lD(a, 0);

  sampleD_kernel2(a);

  cfg->timeD = ((num_t) clock() - start) / (SECONDS * CLOCKS_PER_SEC);
}
