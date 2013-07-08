#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

__host__ Chain *allocChain(Config *cfg){
  Chain *a;
  cudaMalloc((void**) &a, sizeof(Chain));
  allocChain_kernel<<<1, 1>>>(a, cfg->M, cfg->N, cfg->G);
  return a;
}

__global__ void allocChain_kernel(Chain *a, int M, int N, int G){
  int m, n;
  
  /* data */  
  
  a->y = (count_t**) malloc(N * sizeof(count_t*));
  for(n = 0; n < N; ++n)
    a->y[n] = (count_t*) malloc(G * sizeof(count_t));

  a->yMeanG = (num_t*) malloc(N * sizeof(num_t));
  a->grp = (int*) malloc(N * sizeof(int));

  /* parameters */

  a->c      = (num_t**)  malloc((M + 1) * sizeof(num_t*));
  a->sigC   = (num_t*)   malloc((M + 1) * sizeof(num_t));
  a->eps    = (num_t***) malloc((M + 1) * sizeof(num_t**));
  a->eta    = (num_t**)  malloc((M + 1) * sizeof(num_t*));
  a->d      = (num_t*)   malloc((M + 1) * sizeof(num_t));  
  a->tau    = (num_t*)   malloc((M + 1) * sizeof(num_t));
  a->phi    = (num_t**)  malloc((M + 1) * sizeof(num_t*));
  a->thePhi = (num_t*)   malloc((M + 1) * sizeof(num_t));
  a->sigPhi = (num_t*)   malloc((M + 1) * sizeof(num_t));
  a->alp    = (num_t**)  malloc((M + 1) * sizeof(num_t*));
  a->theAlp = (num_t*)   malloc((M + 1) * sizeof(num_t));
  a->sigAlp = (num_t*)   malloc((M + 1) * sizeof(num_t));
  a->piAlp  = (num_t*)   malloc((M + 1) * sizeof(num_t));
  a->del    = (num_t**)  malloc((M + 1) * sizeof(num_t*));
  a->theDel = (num_t*)   malloc((M + 1) * sizeof(num_t));
  a->sigDel = (num_t*)   malloc((M + 1) * sizeof(num_t));
  a->piDel  = (num_t*)   malloc((M + 1) * sizeof(num_t));
  
  for(m = 0; m <= M; ++m){
  
    a->c[m]   = (num_t*)  malloc(N * sizeof(num_t));
    a->eps[m] = (num_t**) malloc(N * sizeof(num_t*)); 
    a->eta[m] = (num_t*)  malloc(G * sizeof(num_t));
    a->phi[m] = (num_t*)  malloc(G * sizeof(num_t));
    a->alp[m] = (num_t*)  malloc(G * sizeof(num_t));
    a->del[m] = (num_t*)  malloc(G * sizeof(num_t));
    
    for(n = 0; n < N; ++n)
      a->eps[m][n] = (num_t*) malloc(G * sizeof(num_t));
  }
  
  /* temporary and return values */
  
  a->tmp1 = (num_t*) malloc(G * sizeof(num_t));
  a->tmp2 = (num_t*) malloc(G * sizeof(num_t));

  a->Old   = (num_t*) malloc(N * sizeof(num_t));
  a->New   = (num_t*) malloc(N * sizeof(num_t));
  a->lOld  = (num_t*) malloc(N * sizeof(num_t));
  a->lNew  = (num_t*) malloc(N * sizeof(num_t));

  /* tuning parameters for Metropolis steps */
  
  a->tuneC = (num_t*) malloc(N * sizeof(num_t));
  a->tunePhi = (num_t*) malloc(G * sizeof(num_t));
  
  a->tuneEps = (num_t**) malloc(N * sizeof(num_t*));
  for(n = 0; n < N; ++n)
    a->tuneEps[n] = (num_t*) malloc(G * sizeof(num_t));

  /* number of acceptances for Metropolis steps */

  a->accC   = (int*) malloc(N * sizeof(int));
  a->accPhi = (int*) malloc(G * sizeof(int));
  a->accAlp = (int*) malloc(G * sizeof(int));
  a->accDel = (int*) malloc(G * sizeof(int));
  
  a->accEps = (int**) malloc(N * sizeof(int*));
  for(n = 0; n < N; ++n)
    a->accEps[n] = (int*) malloc(G * sizeof(int));
}