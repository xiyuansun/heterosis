#include <Chain.h>
#include <Config.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

__host__ Chain *allocChain(Config *cfg, int onHost){

  Chain *a = NULL, *host_a = (Chain*) malloc(sizeof(Chain));
  
  /* data */  
    
  ALLOC(host_a->y, (count_t*), cfg->N * cfg->G * sizeof(count_t), onHost);
  ALLOC(host_a->yMeanG, cfg->N * sizeof(num_t), (num_t*), onHost);
  ALLOC(host_a->grp, cfg->N * sizeof(int), (int*), onHost);

  /* parameters */
  
  ALLOC(host_a->c, (num_t*), (cfg->M + 1) * cfg->N * sizeof(num_t), onHost);
  ALLOC(host_a->sigC, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(host_a->eps, (num_t*), (cfg->M + 1) * cfg->N * cfg->G * sizeof(num_t), onHost);
  ALLOC(host_a->eta, (num_t*), (cfg->M + 1) * cfg->G * sizeof(num_t), onHost);
  ALLOC(host_a->d, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(host_a->tau, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(host_a->phi, (num_t*), (cfg->M + 1) * cfg->G * sizeof(num_t), onHost);
  ALLOC(host_a->thePhi, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(host_a->sigPhi, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(host_a->alp, (num_t*), (cfg->M + 1) * cfg->G * sizeof(num_t), onHost);
  ALLOC(host_a->theAlp, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(host_a->sigAlp, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(host_a->piAlp, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(host_a->del, (num_t*), (cfg->M + 1) * cfg->G * sizeof(num_t), onHost);
  ALLOC(host_a->theDel, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(host_a->sigDel, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(host_a->piDel, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  
  /* temporary and return values */
    
  ALLOC(host_a->tmp1, (num_t*), cfg->G * sizeof(num_t), onHost);
  ALLOC(host_a->tmp2, (num_t*), cfg->G * sizeof(num_t), onHost);
  
  ALLOC(host_a->Old, (num_t*), cfg->N * sizeof(num_t), onHost);
  ALLOC(host_a->New, (num_t*), cfg->N * sizeof(num_t), onHost);
  ALLOC(host_a->lOld, (num_t*), cfg->N * sizeof(num_t), onHost);
  ALLOC(host_a->lNew, (num_t*), cfg->N * sizeof(num_t), onHost);

  /* tuning parameters for Metropolis steps */
  
  ALLOC(host_a->tuneC, (num_t*), cfg->N * sizeof(num_t), onHost);
  ALLOC(host_a->tunePhi, (num_t*), cfg->G * sizeof(num_t), onHost);
  ALLOC(host_a->tuneEps, (num_t*), cfg->N * cfg->G * sizeof(num_t), onHost);

  /* number of acceptances for Metropolis steps */
  
  ALLOC(host_a->accC, (int*), cfg->N * sizeof(int), onHost);
  ALLOC(host_a->accPhi, (int*), cfg->G * sizeof(int), onHost);
  ALLOC(host_a->accAlp, (int*), cfg->G * sizeof(int), onHost);
  ALLOC(host_a->accDel, (int*), cfg->G * sizeof(int), onHost);
  ALLOC(host_a->accEps, (int*), cfg->N * cfg->G * sizeof(int), onHost);
   
  if(onHost){
    a = host_a;
  } else {
    CUDA_CALL(cudaMalloc((void **) &a, sizeof(Chain)));
    CUDA_CALL(cudaMemcpy(a, host_a, sizeof(Chain), cudaMemcpyHostToDevice));
    free(host_a);
  }
    
  return a;
}