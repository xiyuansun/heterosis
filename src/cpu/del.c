#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

num_t delProp(Chain *a, int g){ /* device */    

  num_t gam = a->gamDel;
  num_t sig = a->sigDel;

  num_t gprec = 1/(gam * gam);
  num_t sprec = 1/(sig * sig);

  num_t avg = (a->del[g] * sprec) / (gprec + sprec);
  num_t s = gam * gam + sig * sig;
  num_t u = runiform(0, 1);
  num_t nw;

  if(u < a->piDel){
    nw = 0;
  } else {
    nw = rnormal(avg, s);
  }

  return nw;
}

num_t lDel(Chain *a, int g, num_t arg){ /* device */ 
  int n, G = a->G;
  num_t s = 0, tmp; 
  
  for(n = 0; n < a->N; ++n){
    if(a->grp[n] != 2){
      tmp = mu(a, n, a->phi[g], a->alp[g], arg);
      s += a->y[iG(n, g)] * tmp - exp(a->c[n] + 
          a->eps[iG(n, g)] + tmp);
    }
  }
 
  if(arg * arg > 1e-6){
    tmp = -pow(arg - a->theDel, 2) / (2 * pow(a->sigDel, 2)) -
                log(1 - a->piDel);
  } else {
    tmp = log(a->piDel);
  }

  return s + tmp;
}

void sampleDel_kernel(Chain *a){ /* kernel <<<G, 1>>> */
  int g;
  num_t old, nw, dl, lp, lu;

  for(g = 0; g < a->G; ++g){ 

    old = a->del[g];
    
    if(!a->delPrior){
      nw = delProp(a, g);
      dl = lDel(a, g, nw) - lDel(a, g, old);
    }
    
    lp = 0 < dl? 0 : dl;
    lu = log(runiform(0, 1));
    
    if(lu < lp){ /* accept */
      a->del[g] = nw;
      
      if(a->m > a->burnin)
        ++a->accDel[g];
    } else { /* reject */
      a->del[g] = old;
    }
  }
}

void sampleDel(Chain *a, Config *cfg){ /* host */

  clock_t start = clock();
  
  if(cfg->verbose)
    printf("del ");
  
  sampleDel_kernel(a);

  cfg->timeDel = ((num_t) clock() - start) / (SECONDS * CLOCKS_PER_SEC);
}
