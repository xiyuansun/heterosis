#include <Chain.h>
#include <Config.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

__host__ void freeChain(Chain *a, Config *cfg, int onHost){

  if(cfg->verbose && !onHost)
    printf("  Freeing chain.\n\n");

  FREE(a->states, onHost);
  
  /* data */  

  FREE(a->y, onHost);
  FREE(a->yMeanG, onHost);
  FREE(a->grp, onHost);

  /* parameters */
  
  FREE(a->c, onHost);
  FREE(a->eps, onHost);
  FREE(a->eta, onHost);
  FREE(a->phi, onHost);
  FREE(a->alp, onHost);
  FREE(a->del, onHost);

  /* temporary and return values */
  
  FREE(a->tmp1, onHost);
  FREE(a->tmp2, onHost);

  FREE(a->Old, onHost);
  FREE(a->New, onHost);
  FREE(a->lOld, onHost);
  FREE(a->lNew, onHost);
  
  /* tuning parameters for Metropolis steps */
  
  FREE(a->tuneC, onHost);
  FREE(a->tunePhi, onHost);
  FREE(a->tuneEps, onHost);

  /* number of acceptances for Metropolis steps */

  FREE(a->accC, onHost);
  FREE(a->accPhi, onHost);
  FREE(a->accAlp, onHost);
  FREE(a->accDel, onHost);
  FREE(a->accEps, onHost);
  
  /* counters towards heterosis and differential expression */

  FREE(a->dex, onHost);
  FREE(a->hph, onHost);
  FREE(a->lph, onHost);
  FREE(a->mph, onHost);
  
  /* for computing DIC */
  
  FREE(a->meanC, onHost);
  FREE(a->meanPhi, onHost);  
  FREE(a->meanAlp, onHost);
  FREE(a->meanDel, onHost);
  FREE(a->meanEps, onHost);  
  
  free(a);
}
