#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

__host__ Chain *allocChainHost(Config *cfg){
  int N = cfg->N, G = cfg->G;
  Chain *a = (Chain*) malloc(sizeof(Chain));
  
  /* curand states */
  
  a->states = (curandState_t*) malloc(MAX_NG * sizeof(curandState_t));

  /* data */   
  
  a->y = (count_t*) malloc(cfg->N * cfg->G * sizeof(count_t));
  a->yMeanG = (num_t*) malloc(cfg->N * sizeof(num_t));
  a->grp = (int*) malloc(cfg->N * sizeof(int));

  /* parameters */

  a->c      = (num_t*) malloc( cfg->N * sizeof(num_t));
  a->eps    = (num_t*) malloc( cfg->N * cfg->G * sizeof(num_t));
  a->eta    = (num_t*) malloc( cfg->G * sizeof(num_t));
  a->phi    = (num_t*) malloc( cfg->G * sizeof(num_t));
  a->alp    = (num_t*) malloc( cfg->G * sizeof(num_t));
  a->del    = (num_t*) malloc( cfg->G * sizeof(num_t));
  
  /* temporary and return values */

  a->tmp1 = (num_t*) malloc(cfg->G * sizeof(num_t));
  a->tmp2 = (num_t*) malloc(cfg->G * sizeof(num_t));

  a->Old   = (num_t*) malloc(cfg->N * sizeof(num_t));
  a->New   = (num_t*) malloc(cfg->N * sizeof(num_t));
  a->lOld  = (num_t*) malloc(cfg->N * sizeof(num_t));
  a->lNew  = (num_t*) malloc(cfg->N * sizeof(num_t));

  /* tuning parameters for Metropolis steps */
  
  a->tuneC = (num_t*) malloc(cfg->N * sizeof(num_t));
  a->tunePhi = (num_t*) malloc(cfg->G * sizeof(num_t));
  a->tuneEps = (num_t*) malloc(cfg->N * cfg->G * sizeof(num_t));

  /* number of acceptances for Metropolis steps */

  a->accC   = (int*) malloc(cfg->N * sizeof(int));
  a->accPhi = (int*) malloc(cfg->G * sizeof(int));
  a->accAlp = (int*) malloc(cfg->G * sizeof(int));
  a->accDel = (int*) malloc(cfg->G * sizeof(int));
  a->accEps = (int*) malloc(cfg->N * cfg->G * sizeof(int));
  
  /* counts toward differential expression and heterosis */
  
  a->dex = (int*) malloc(cfg->G * sizeof(int));
  a->hph = (int*) malloc(cfg->G * sizeof(int));
  a->lph = (int*) malloc(cfg->G * sizeof(int));
  a->mph = (int*) malloc(cfg->G * sizeof(int));
    
  return a;
}

#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

__host__ void allocChainDevice(Chain **host_a, Chain **dev_a, Config *cfg){

  *host_a = (Chain*) malloc(sizeof(Chain));

  /* curand states */
  
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->states), cfg->N * cfg->G * sizeof(curandState_t)));
  
  /* data */   
  
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->y), cfg->N * cfg->G * sizeof(count_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->yMeanG), cfg->N * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->grp), cfg->N * sizeof(int)));

  /* parameters */

  CUDA_CALL(cudaMalloc((void**) &((*host_a)->c), cfg->N * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->eps), cfg->N * cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->eta), cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->phi), cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->alp), cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->del), cfg->G * sizeof(num_t)));
  
  /* temporary and return values */

  CUDA_CALL(cudaMalloc((void**) &((*host_a)->tmp1), cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->tmp2), cfg->G * sizeof(num_t)));

  CUDA_CALL(cudaMalloc((void**) &((*host_a)->Old), cfg->N * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->New), cfg->N * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->lOld), cfg->N * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->lNew), cfg->N * sizeof(num_t)));

  /* tuning parameters for Metropolis steps */
  
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->tuneC), cfg->N * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->tunePhi), cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->tuneEps), cfg->N * cfg->G * sizeof(num_t)));

  /* number of acceptances for Metropolis steps */

  CUDA_CALL(cudaMalloc((void**) &((*host_a)->accC), cfg->N * sizeof(int)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->accPhi), cfg->G * sizeof(int)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->accAlp),cfg->G * sizeof(int)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->accDel), cfg->G * sizeof(int)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->accEps), cfg->N * cfg->G * sizeof(int)));
  
  /* counts toward differential expression and heterosis */
  
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->dex), cfg->G * sizeof(int)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->hph), cfg->G * sizeof(int)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->lph), cfg->G * sizeof(int)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->mph), cfg->G * sizeof(int)));

  /* for computing DIC */
  
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->meanC), cfg->N * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->meanPhi), cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->meanAlp), cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->meanDel), cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->meanEps), cfg->N * cfg->G * sizeof(num_t)));  
  
  /* pointer to chain on the device */

  CUDA_CALL(cudaMalloc((void**) dev_a, sizeof(Chain)));
  CUDA_CALL(cudaMemcpy(*dev_a, *host_a, sizeof(Chain), cudaMemcpyHostToDevice));
}