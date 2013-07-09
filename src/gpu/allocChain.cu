#include <Chain.h>
#include <Config.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

__host__ void allocChain(Chain *a, Config *cfg){

  Chain *host_a = (Chain*) malloc(sizeof(Chain));
  a = host_a;
  CUDA_CALL(cudaMalloc((void **) &a, sizeof(Chain)));
  
  /* data */  
  
  CUDA_CALL(cudaMalloc((void **) &(host_a->y), cfg->N * cfg->G * sizeof(count_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->yMeanG), cfg->N * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->grp), cfg->N * sizeof(int)));

  /* parameters */

  CUDA_CALL(cudaMalloc((void **) &(host_a->c), (cfg->M + 1) * cfg->N * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->sigC), (cfg->M + 1) * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->eps), (cfg->M + 1) * cfg->N * cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->eta), (cfg->M + 1) * cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->d), (cfg->M + 1) * sizeof(num_t)));  
  CUDA_CALL(cudaMalloc((void **) &(host_a->tau), (cfg->M + 1) * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->phi), (cfg->M + 1) * cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->thePhi), (cfg->M + 1) * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->sigPhi), (cfg->M + 1) * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->alp), (cfg->M + 1) * cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->theAlp), (cfg->M + 1) * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->sigAlp), (cfg->M + 1) * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->piAlp), (cfg->M + 1) * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->del), (cfg->M + 1) * cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->theDel), (cfg->M + 1) * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->sigDel), (cfg->M + 1) * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->piDel), (cfg->M + 1) * sizeof(num_t)));
  
  /* temporary and return values */
  
  CUDA_CALL(cudaMalloc((void **) &(host_a->tmp1), cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->tmp2), cfg->G * sizeof(num_t)));

  CUDA_CALL(cudaMalloc((void **) &(host_a->Old), cfg->N * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->New), cfg->N * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->lOld), cfg->N * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->lNew), cfg->N * sizeof(num_t)));

  /* tuning parameters for Metropolis steps */
  
  CUDA_CALL(cudaMalloc((void **) &(host_a->tuneC), cfg->N * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->tunePhi), cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->tuneEps), cfg->N * cfg->G * sizeof(num_t)));

  /* number of acceptances for Metropolis steps */

  CUDA_CALL(cudaMalloc((void **) &(host_a->accC), cfg->N * sizeof(int)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->accPhi), cfg->G * sizeof(int)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->accAlp), cfg->G * sizeof(int)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->accDel), cfg->G * sizeof(int)));
  CUDA_CALL(cudaMalloc((void **) &(host_a->accEps), cfg->N * cfg->G * sizeof(int)));
    
  CUDA_CALL(cudaMemcpy(a, host_a, sizeof(Chain), cudaMemcpyHostToDevice));
  free(host_a);
}