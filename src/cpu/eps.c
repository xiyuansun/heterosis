#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

num_t lEps(Chain *a, int n, int g, num_t arg){ /* device */
  int N = a->N, G = a->G;
  return a->y[iG(n, g)] * arg - exp(a->c[n] + arg + mu(a, n, a->phi[g], 
                                     a->alp[g], a->del[g])) 
                          - (arg * arg) / (2 * pow(a->eta[g], 2));
}

void sampleEps_kernel(Chain *a){ /* kernel <<<N, G>>> */
  int n, g, N = a->N, G = a->G;
  num_t old, nw, dl, lp, lu;

  for(g = 0; g < a->G; ++g){
    for(n = 0; n < a->N; ++n){ 
      old = a->eps[iG(n, g)];
      nw = rnormal(old, a->tuneEps[iG(n, g)]);

      dl = lEps(a, n, g, nw) - lEps(a, n, g, old);
      lp = 0 < dl ? 0 : dl;
      lu = log(runiform(0, 1));
      
      if(lu < lp){ /* accept */
        a->eps[iG(n, g)] = nw;
        a->tuneEps[iG(n, g)] *= 1.1;
        
        if(a->m >= a->burnin)
          ++a->accEps[iG(n, g)]; 
      } else { /* reject */
        a->eps[iG(n, g)] = old;
        a->tuneEps[iG(n, g)] /= 1.1;
      }
    }
  }
}

void sampleEps(Chain *a, Config *cfg){ /* host */

  clock_t start = clock();

  if(cfg->verbose)
    printf("eps ");
  
  sampleEps_kernel(a);

  cfg->timeEps = ((num_t) clock() - start) / (SECONDS * CLOCKS_PER_SEC);
}