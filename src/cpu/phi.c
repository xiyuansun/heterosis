#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

num_t lPhi(Chain *a, int g, num_t arg){ /* device */
  int n, G = a->G;
  num_t ret, s = 0, tmp = 0; 

  for(n = 0; n < a->N; ++n){
    tmp = mu(a, n, arg, a->alp[g], a->del[g]);
    s += a->y[iG(n, g)] * tmp - exp(a->c[n] + 
        a->eps[iG(n, g)] + tmp);
  }
 
  ret = s - pow(arg - a->thePhi, 2) / (2 * pow(a->sigPhi, 2));
  return ret;
}

void samplePhi_kernel(Chain *a){ /* kernel <<<G, 1>>> */
  int g;
  num_t old, nw, dl, lp, lu;
  
  for(g = 0; g < a->G; ++g){ 

    old = a->phi[g];
    nw = rnormal(old, a->tunePhi[g]);

    dl = lPhi(a, g, nw) - lPhi(a, g, old);
    lp = 0 < dl ? 0 : dl;
    lu = log(runiform(0, 1)); 
    
    if(lu < lp){ /* accept */
      a->phi[g] = nw;
      a->tunePhi[g] *= 1.1; 
      
      if(a->m > a->burnin)
        ++a->accPhi[g];
    } else { /* reject */
      a->phi[g] = old;
      a->tunePhi[g] /= 1.1; 
    }
    
    if(a->m > a->burnin)
      a->sumPhi[g] += a->phi[g];
  }
}

void samplePhi(Chain *a, Config *cfg){ /* host */

  clock_t start = clock();
  
  if(cfg->verbose)
    printf("phi ");
  
  samplePhi_kernel(a);

  cfg->timePhi = ((num_t) clock() - start) / (SECONDS * CLOCKS_PER_SEC);
}
