#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

num_t alpProp(Chain *a, int g){ /* device */

  num_t gam = a->gamAlp;
  num_t sig = a->sigAlp;

  num_t gprec = 1/(gam * gam);
  num_t sprec = 1/(sig * sig);

  num_t avg = (a->alp[g] * sprec) / (gprec + sprec);
  num_t s = gam * gam + sig * sig;
  num_t u = runiform(0, 1);
  num_t nw;
  
  if(u < a->piAlp){
    nw = 0;
  } else {
    nw = rnormal(avg, s);
  }

  return nw;
}

num_t lAlp(Chain *a, int g, num_t arg){ /* device */
  
  int n, G = a->G;
  num_t s = 0, tmp;
   
  for(n = 0; n < a->N; ++n){
    if(a->grp[n] != 2){
      tmp = mu(a, n, a->phi[g], arg, a->del[g]);
      s += a->y[iG(n, g)] * tmp - exp(a->c[n] + 
          a->eps[iG(n, g)] + tmp);
    }
  }
 
  if(arg * arg > 1e-6){
    tmp = -pow(arg - a->theAlp, 2) / (2 * pow(a->sigAlp, 2)) -
                log(1 - a->piAlp);
  } else {
    tmp = log(a->piAlp);
  }

  return s + tmp;
}

void sampleAlp_kernel(Chain *a){ /* kernel <<<G, 1>>> */

  int g;
  num_t old, nw, dl, lp, lu;

  for(g = 0; g < a->G; ++g){ 

    old = a->alp[g];
    
    if(!a->alpPrior){
      nw = alpProp(a, g);
      dl = lAlp(a, g, nw) - lAlp(a, g, old);
    }

    lp = 0 < dl ? 0 : dl;
    lu = log(runiform(0, 1));
    
    if(lu < lp){ /* accept */
      a->alp[g] = nw;
      
      if(a->m > a->burnin)
        ++a->accAlp[g];
    } 
  }
}

void sampleAlp(Chain *a, Config *cfg){ /* host */

  clock_t start = clock();
  
  if(cfg->verbose)
    printf("alpha ");
    
  sampleAlp_kernel(a);

  cfg->timeAlp = ((num_t) clock() - start) / (SECONDS * CLOCKS_PER_SEC);
}