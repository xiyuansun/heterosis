#include <Chain.h>
#include <Config.h>
#include <functions.h>
#include <numericTypes.h>
#include <stdio.h>
#include <stdlib.h>

void freeChainHost(Chain *a, Config *cfg){
  
  int m, n;
  
  /* data */  
  
  for(n = 0; n < cfg->N; ++n)
    free(a->y[n]);

  free(a->y);
  free(a->yMeanG);
  free(a->grp);

  /* parameters */
  
  for(m = 0; m <= cfg->M; ++m){

    for(n = 0; n < cfg->N; ++n)
      free(a->eps[m][n]);
      
    free(a->c[m]);
    free(a->eps[m]);
    free(a->eta[m]);
    free(a->phi[m]);
    free(a->alp[m]);
    free(a->del[m]);
  }
  
  free(a->c);
  free(a->sigC);
  free(a->eps);
  free(a->eta);
  free(a->d);
  free(a->tau);
  free(a->phi);
  free(a->thePhi);
  free(a->sigPhi);
  free(a->alp);
  free(a->theAlp);
  free(a->sigAlp);
  free(a->piAlp);
  free(a->del);
  free(a->theDel);
  free(a->sigDel);
  free(a->piDel);
  
  /* temporary and return values */
  
  for(n = 0; n < cfg->N; ++n){
    free(a->Old[n]);
    free(a->New[n]);
    free(a->lOld[n]);
    free(a->lNew[n]);
  }
  
  free(a->tmp1);
  free(a->tmp2);

  free(a->Old);
  free(a->New);
  free(a->lOld);
  free(a->lNew);
  
  /* tuning parameters for Metropolis steps */

  for(n = 0; n < cfg->N; ++n)
    free(a->tuneEps[n]);
  
  free(a->tuneC);
  free(a->tunePhi);
  free(a->tuneEps);

  /* number of acceptances for Metropolis steps */

  for(n = 0; n < cfg->N; ++n)
    free(a->accEps[n]);

  free(a->accC);
  free(a->accPhi);
  free(a->accAlp);
  free(a->accDel);
  free(a->accEps);
}