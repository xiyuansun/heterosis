#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

Chain *allocChain(Config *cfg){
  int m, n;
  Chain *a;
  
  a = (Chain*) malloc(sizeof(Chain));
  
  /* data */  
  
  a->y = (count_t**) malloc(cfg->N * sizeof(count_t*));
  for(n = 0; n < cfg->N; ++n)
    a->y[n] = (count_t*) malloc(cfg->G * sizeof(count_t));

  a->yMeanG = (num_t*) malloc(cfg->N * sizeof(num_t));
  a->grp = (int*) malloc(cfg->N * sizeof(int));

  /* parameters */

  a->c      = (num_t**)  malloc((cfg->M + 1) * sizeof(num_t*));
  a->sigC   = (num_t*)   malloc((cfg->M + 1) * sizeof(num_t));
  a->eps    = (num_t***) malloc((cfg->M + 1) * sizeof(num_t**));
  a->eta    = (num_t**)  malloc((cfg->M + 1) * sizeof(num_t*));
  a->d      = (num_t*)   malloc((cfg->M + 1) * sizeof(num_t));  
  a->tau    = (num_t*)   malloc((cfg->M + 1) * sizeof(num_t));
  a->phi    = (num_t**)  malloc((cfg->M + 1) * sizeof(num_t*));
  a->thePhi = (num_t*)   malloc((cfg->M + 1) * sizeof(num_t));
  a->sigPhi = (num_t*)   malloc((cfg->M + 1) * sizeof(num_t));
  a->alp    = (num_t**)  malloc((cfg->M + 1) * sizeof(num_t*));
  a->theAlp = (num_t*)   malloc((cfg->M + 1) * sizeof(num_t));
  a->sigAlp = (num_t*)   malloc((cfg->M + 1) * sizeof(num_t));
  a->piAlp  = (num_t*)   malloc((cfg->M + 1) * sizeof(num_t));
  a->del    = (num_t**)  malloc((cfg->M + 1) * sizeof(num_t*));
  a->theDel = (num_t*)   malloc((cfg->M + 1) * sizeof(num_t));
  a->sigDel = (num_t*)   malloc((cfg->M + 1) * sizeof(num_t));
  a->piDel  = (num_t*)   malloc((cfg->M + 1) * sizeof(num_t));
  
  for(m = 0; m <= cfg->M; ++m){
  
    a->c[m]   = (num_t*)  malloc(cfg->N * sizeof(num_t));
    a->eps[m] = (num_t**) malloc(cfg->N * sizeof(num_t*)); 
    a->eta[m] = (num_t*)  malloc(cfg->G * sizeof(num_t));
    a->phi[m] = (num_t*)  malloc(cfg->G * sizeof(num_t));
    a->alp[m] = (num_t*)  malloc(cfg->G * sizeof(num_t));
    a->del[m] = (num_t*)  malloc(cfg->G * sizeof(num_t));
    
    for(n = 0; n < cfg->N; ++n)
      a->eps[m][n] = (num_t*) malloc(cfg->G * sizeof(num_t));
  }
  
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
  
  a->tuneEps = (num_t**) malloc(cfg->N * sizeof(num_t*));
  for(n = 0; n < cfg->N; ++n)
    a->tuneEps[n] = (num_t*) malloc(cfg->G * sizeof(num_t));

  /* number of acceptances for Metropolis steps */

  a->accC   = (int*) malloc(cfg->N * sizeof(int));
  a->accPhi = (int*) malloc(cfg->G * sizeof(int));
  a->accAlp = (int*) malloc(cfg->G * sizeof(int));
  a->accDel = (int*) malloc(cfg->G * sizeof(int));
  
  a->accEps = (int**) malloc(cfg->N * sizeof(int*));
  for(n = 0; n < cfg->N; ++n)
    a->accEps[n] = (int*) malloc(cfg->G * sizeof(int));
    
  return a;
}