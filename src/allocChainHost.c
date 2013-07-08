#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

Chain *allocChainHost(Config *cfg){
  int m, n;
  Chain *a;
  
  a = malloc(sizeof(Chain));
  
  /* data */  
  
  a->y = malloc(cfg->N * sizeof(count_t*));
  for(n = 0; n < cfg->N; ++n)
    a->y[n] = malloc(cfg->G * sizeof(count_t));

  a->yMeanG = malloc(cfg->N * sizeof(count_t));
  a->grp = malloc(cfg->N * sizeof(count_t));

  /* parameters */

  a->c      = malloc((cfg->M + 1) * sizeof(num_t*));
  a->sigC   = malloc((cfg->M + 1) * sizeof(num_t));
  a->eps    = malloc((cfg->M + 1) * sizeof(num_t**));
  a->eta    = malloc((cfg->M + 1) * sizeof(num_t*));
  a->d      = malloc((cfg->M + 1) * sizeof(num_t));  
  a->tau    = malloc((cfg->M + 1) * sizeof(num_t));
  a->phi    = malloc((cfg->M + 1) * sizeof(num_t*));
  a->thePhi = malloc((cfg->M + 1) * sizeof(num_t));
  a->sigPhi = malloc((cfg->M + 1) * sizeof(num_t));
  a->alp    = malloc((cfg->M + 1) * sizeof(num_t*));
  a->theAlp = malloc((cfg->M + 1) * sizeof(num_t));
  a->sigAlp = malloc((cfg->M + 1) * sizeof(num_t));
  a->piAlp  = malloc((cfg->M + 1) * sizeof(num_t));
  a->del    = malloc((cfg->M + 1) * sizeof(num_t*));
  a->theDel = malloc((cfg->M + 1) * sizeof(num_t));
  a->sigDel = malloc((cfg->M + 1) * sizeof(num_t));
  a->piDel  = malloc((cfg->M + 1) * sizeof(num_t));
  
  for(m = 0; m <= cfg->M; ++m){
  
    a->c[m]   = malloc(cfg->N * sizeof(num_t));
    a->eps[m] = malloc(cfg->N * sizeof(num_t*)); 
    a->eta[m] = malloc(cfg->G * sizeof(num_t));
    a->phi[m] = malloc(cfg->G * sizeof(num_t));
    a->alp[m] = malloc(cfg->G * sizeof(num_t));
    a->del[m] = malloc(cfg->G * sizeof(num_t));
    
    for(n = 0; n < cfg->N; ++n)
      a->eps[m][n] = malloc(cfg->G * sizeof(num_t));
  }
  
  /* temporary and return values */
  
  a->tmp1 = malloc(mx * sizeof(num_t));
  a->tmp2 = malloc(mx * sizeof(num_t));

  a->Old   = malloc(cfg->N * sizeof(num_t));
  a->New   = malloc(cfg->N * sizeof(num_t));
  a->lOld  = malloc(cfg->N * sizeof(num_t));
  a->lNew  = malloc(cfg->N * sizeof(num_t));

  /* tuning parameters for Metropolis steps */
  
  a->tuneC = malloc(cfg->N * sizeof(num_t));
  a->tunePhi = malloc(cfg->G * sizeof(num_t));
  
  a->tuneEps = malloc(cfg->N * sizeof(num_t*));
  for(n = 0; n < cfg->N; ++n)
    a->tuneEps[n] = malloc(cfg->G * sizeof(num_t));

  /* number of acceptances for Metropolis steps */

  a->accC   = malloc(cfg->N * sizeof(num_t));
  a->accPhi = malloc(cfg->G * sizeof(num_t));
  a->accAlp = malloc(cfg->G * sizeof(num_t));
  a->accDel = malloc(cfg->G * sizeof(num_t));
  
  a->accEps = malloc(cfg->N * sizeof(num_t*));
  for(n = 0; n < cfg->N; ++n)
    a->accEps[n] = malloc(cfg->G * sizeof(num_t));
    
  return a;
}