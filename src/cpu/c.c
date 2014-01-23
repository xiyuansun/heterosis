#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void lC_kernel1(Chain *a, int n){ /* kernel <<<G, 1>>> */
  int g, G = a->G;
  
  for(g = 0; g < a->G; ++g)
    a->tmp1[g] = exp(a->eps[iG(n, g)] + mu(a, n, a->phi[g], 
                    a->alp[g], a->del[g]));
}

void lC_kernel2(Chain *a, int n){ /* parallel pairwise sum */
  int g;
  a->tmp2[0] = 0;
  
  for(g = 0; g < a->G; ++g)
    a->tmp2[0] += a->tmp1[g];
}

void lC_kernel3(Chain *a, int n, int newArg){ /* kernel <<<1, 1>>> */
  num_t arg, ret;

  if(newArg){
    arg = a->New[n];
  } else {
    arg = a->Old[n];
  }

  ret = arg * a->G * a->yMeanG[n] - exp(arg) * a->tmp2[0] - (arg*arg) / 
        (2 * a->sigC * a->sigC);

  if(newArg){
    a->lNew[n] = ret;
  } else {
    a->lOld[n] = ret;
  }
}

void lC(Chain *a, int n, int newArg){ /* host */
  lC_kernel1(a, n);
  lC_kernel2(a, n);
  lC_kernel3(a, n, newArg);
}

void sampleC_kernel1(Chain *a){ /* kernel <<<1, N>>> */
  int n;
  
  for(n = 0; n < a->N; ++n){
    a->Old[n] = a->c[n];
    a->New[n] = rnormal(a->Old[n], a->tuneC[n]);
    
  }
}

void sampleC_kernel2(Chain *a){ /* kernel <<<1, N>>> */
  int n;
  num_t dl, lp, lu;

  for(n = 0; n < a->N; ++n){ 

    dl = a->lNew[n] - a->lOld[n];
    lp = 0 < dl ? 0 : dl;
    lu = log(runiform(0, 1));
      
    if(lu < lp){ /* accept */
      a->c[n] = a->New[n];
      a->tuneC[n] *= 1.1; /* Increase the proposal variance to avoid  
                                       gettiG stuck in a mode */
                                       
      if(a->m > a->burnin)                                 
        ++a->accC[n];
    } else { /* reject */
      a->tuneC[n] /= 1.1; /* If you're rejectiG too often, decrease the  
                                       proposal variance to sample closer to 
                                       the last accepted value. */
    }
  }
}

void sampleC(Chain *a, Config *cfg){ /* host */
  int n;
  clock_t start = clock();

  if(cfg->verbose)
    printf("c ");
  
  sampleC_kernel1(a);

  for(n = 0; n < a->N; ++n){ 
    lC(a, n, 1);
    lC(a, n, 0);
  }

  sampleC_kernel2(a);
    
  cfg->timeC = ((num_t) clock() - start) / (SECONDS * CLOCKS_PER_SEC);
}
