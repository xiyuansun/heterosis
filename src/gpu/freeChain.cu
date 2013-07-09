#include <Chain.h>
#include <Config.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

__host__ void freeChain(Chain *a, Config *cfg){
  
  Chain *host_a = (Chain*) malloc(sizeof(Chain));
  CUDA_CALL(cudaMemcpy(host_a, a, sizeof(Chain), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaFree(a));
  
  /* data */  

  CUDA_CALL(cudaFree(host_a->y));
  CUDA_CALL(cudaFree(host_a->yMeanG));
  CUDA_CALL(cudaFree(host_a->grp));

  /* parameters */
  
  CUDA_CALL(cudaFree(host_a->c));
  CUDA_CALL(cudaFree(host_a->sigC));
  CUDA_CALL(cudaFree(host_a->eps));
  CUDA_CALL(cudaFree(host_a->eta));
  CUDA_CALL(cudaFree(host_a->d));
  CUDA_CALL(cudaFree(host_a->tau));
  CUDA_CALL(cudaFree(host_a->phi));
  CUDA_CALL(cudaFree(host_a->thePhi));
  CUDA_CALL(cudaFree(host_a->sigPhi));
  CUDA_CALL(cudaFree(host_a->alp));
  CUDA_CALL(cudaFree(host_a->theAlp));
  CUDA_CALL(cudaFree(host_a->sigAlp));
  CUDA_CALL(cudaFree(host_a->piAlp));
  CUDA_CALL(cudaFree(host_a->del));
  CUDA_CALL(cudaFree(host_a->theDel));
  CUDA_CALL(cudaFree(host_a->sigDel));
  CUDA_CALL(cudaFree(host_a->piDel));
  
  /* temporary and return values */
  
  CUDA_CALL(cudaFree(host_a->tmp1));
  CUDA_CALL(cudaFree(host_a->tmp2));

  CUDA_CALL(cudaFree(host_a->Old));
  CUDA_CALL(cudaFree(host_a->New));
  CUDA_CALL(cudaFree(host_a->lOld));
  CUDA_CALL(cudaFree(host_a->lNew));
  
  /* tuning parameters for Metropolis steps */
  
  CUDA_CALL(cudaFree(host_a->tuneC));
  CUDA_CALL(cudaFree(host_a->tunePhi));
  CUDA_CALL(cudaFree(host_a->tuneEps));

  /* number of acceptances for Metropolis steps */

  CUDA_CALL(cudaFree(host_a->accC));
  CUDA_CALL(cudaFree(host_a->accPhi));
  CUDA_CALL(cudaFree(host_a->accAlp));
  CUDA_CALL(cudaFree(host_a->accDel));
  CUDA_CALL(cudaFree(host_a->accEps));
  
  free(host_a);
}