#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

num_t lPhiAlpDelJoint(Chain *a, int g, num_t argPhi, num_t argAlp, num_t argDel){ /* device */
 
  int n, G = a->G;
  num_t ret, s = 0, tmp = 0;

  for(n = 0; n < a->N; ++n){
    tmp = mu(a, n, argPhi, argAlp, argDel);
    s += a->y[iG(n, g)] * tmp - exp(a->c[n] + 
         a->eps[iG(n, g)] + tmp);
  }

  /* phi part */
  
  if(!a->phiPrior){
    ret = s - pow(argPhi - a->thePhi, 2) / (2 * pow(a->sigPhi, 2));
  }
  
  /* alpha part */
  
  if(!a->alpPrior){
	if(argAlp * argAlp > 1e-6){
	  tmp = -pow(argAlp - a->theAlp, 2) / (2 * pow(a->sigAlp, 2)) -
				  log(1 - a->piAlp);
	} else {
	  tmp = log(a->piAlp);
	}
  
    ret += tmp;
  }
  
  /* delta part */
  
  if(!a->delPrior){
	if(argDel * argDel > 1e-6){
	  tmp = -pow(argDel - a->theDel, 2) / (2 * pow(a->sigDel, 2)) -
				  log(1 - a->piDel);
	} else {
	  tmp = log(a->piDel);
	}
  
	ret += tmp;
  }
  
  return ret;
}

void samplePhiAlpDelJoint_kernel(Chain *a){ /* kernel <<<G, 1>>> */
  int g;
  num_t oldPhi, newPhi, oldAlp, newAlp, oldDel, newDel;
  num_t dl, lp, lu;

  for(g = 0; g < a->G; ++g){

    oldPhi = a->phi[g];
    newPhi = rnormal(oldPhi, a->tunePhi[g]);

    oldAlp = a->alp[g];
    newAlp = alpProp(a, g);

    oldDel = a->del[g];
    newDel = delProp(a, g);

    dl = lPhiAlpDelJoint(a, g, newPhi, newAlp, newDel) 
       - lPhiAlpDelJoint(a, g, oldPhi, oldAlp, oldDel); 
    lp = 0 < dl ? 0 : dl;
    lu = log(runiform(0, 1));

    if(lu < lp){ /* accept */
      a->phi[g] = newPhi;
      a->alp[g] = newAlp;
      a->del[g] = newDel;

      a->tunePhi[g] *= 1.1; 

      if(a->m > a->burnin){
        ++a->accPhi[g];
        ++a->accAlp[g];
        ++a->accDel[g];
      }
    } else { /* reject */
      a->tunePhi[g] /= 1.1;
    }
  }
}

void samplePhiAlpDelJoint(Chain *a, Config *cfg){ /* host */
  num_t time;
  clock_t start = clock();

  if(cfg->verbose)
    printf("phiAlpDelJoint ");  
  
  samplePhiAlpDelJoint_kernel(a); 
  time = ((num_t) clock() - start) / (SECONDS * CLOCKS_PER_SEC);
  
  cfg->timePhi = time;
  cfg->timeAlp = time;
  cfg->timeDel = time;  
}