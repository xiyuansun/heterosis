#include <Chain.h>
#include <Config.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

__host__ void freeChain(Chain *a, Config *cfg, int onHost){
  
  /* data */  

  FREE(a->y, onHost);
  FREE(a->yMeanG, onHost);
  FREE(a->grp, onHost);

  /* parameters */
  
  FREE(a->c, onHost);
  FREE(a->sigC, onHost);
  FREE(a->eps, onHost);
  FREE(a->eta, onHost);
  FREE(a->d, onHost);
  FREE(a->tau, onHost);
  FREE(a->phi, onHost);
  FREE(a->thePhi, onHost);
  FREE(a->sigPhi, onHost);
  FREE(a->alp, onHost);
  FREE(a->theAlp, onHost);
  FREE(a->sigAlp, onHost);
  FREE(a->piAlp, onHost);
  FREE(a->del, onHost);
  FREE(a->theDel, onHost);
  FREE(a->sigDel, onHost);
  FREE(a->piDel, onHost);
  
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
  
  free(a);
}