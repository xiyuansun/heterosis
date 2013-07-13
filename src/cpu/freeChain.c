#include <Chain.h>
#include <Config.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

void freeChain(Chain *a, Config *cfg){

  if(cfg->verbose)
    printf("  Freeing chain.\n\n");
  
  /* data */  

  free(a->y);
  free(a->yMeanG);
  free(a->grp);

  /* parameters */
  
  free(a->c);
  free(a->eps);
  free(a->eta);
  free(a->phi);
  free(a->alp);
  free(a->del);

  /* temporary and return values */
  
  free(a->tmp1);
  free(a->tmp2);

  free(a->Old);
  free(a->New);
  free(a->lOld);
  free(a->lNew);
  
  /* tuning parameters for Metropolis steps */
  
  free(a->tuneC);
  free(a->tunePhi);
  free(a->tuneEps);

  /* number of acceptances for Metropolis steps */

  free(a->accC);
  free(a->accPhi);
  free(a->accAlp);
  free(a->accDel);
  free(a->accEps);
  
  free(a);
}