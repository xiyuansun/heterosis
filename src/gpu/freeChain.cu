#include <Chain.h>
#include <Config.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

__host__ void freeChain(Chain *a, Config *cfg, int onHost){
  
  Chain *host_a = onHost ? a : (Chain*) malloc(sizeof(Chain));
  
  if(!onHost) {
    CUDA_CALL(cudaMemcpy(host_a, a, sizeof(Chain), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaFree(a)); 
  }
  
  /* data */  

  FREE(host_a->y, onHost);
  FREE(host_a->yMeanG, onHost);
  FREE(host_a->grp, onHost);

  /* parameters */
  
  FREE(host_a->c, onHost);
  FREE(host_a->sigC, onHost);
  FREE(host_a->eps, onHost);
  FREE(host_a->eta, onHost);
  FREE(host_a->d, onHost);
  FREE(host_a->tau, onHost);
  FREE(host_a->phi, onHost);
  FREE(host_a->thePhi, onHost);
  FREE(host_a->sigPhi, onHost);
  FREE(host_a->alp, onHost);
  FREE(host_a->theAlp, onHost);
  FREE(host_a->sigAlp, onHost);
  FREE(host_a->piAlp, onHost);
  FREE(host_a->del, onHost);
  FREE(host_a->theDel, onHost);
  FREE(host_a->sigDel, onHost);
  FREE(host_a->piDel, onHost);
  
  /* temporary and return values */
  
  FREE(host_a->tmp1, onHost);
  FREE(host_a->tmp2, onHost);

  FREE(host_a->Old, onHost);
  FREE(host_a->New, onHost);
  FREE(host_a->lOld, onHost);
  FREE(host_a->lNew, onHost);
  
  /* tuning parameters for Metropolis steps */
  
  FREE(host_a->tuneC, onHost);
  FREE(host_a->tunePhi, onHost);
  FREE(host_a->tuneEps, onHost);

  /* number of acceptances for Metropolis steps */

  FREE(host_a->accC, onHost);
  FREE(host_a->accPhi, onHost);
  FREE(host_a->accAlp, onHost);
  FREE(host_a->accDel, onHost);
  FREE(host_a->accEps, onHost);
  
  free(host_a);
}