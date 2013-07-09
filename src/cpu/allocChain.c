#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

Chain *allocChain(Config *cfg){

  Chain *a = malloc(sizeof(Chain));
  
  /* data */  
  
  a->y = malloc(cfg->N * cfg->G * sizeof(count_t));
  a->yMeanG = malloc(cfg->N * sizeof(num_t));
  a->grp = malloc(cfg->N * sizeof(int));

  /* parameters */

  a->c      = malloc((cfg->M + 1) * cfg->N * sizeof(num_t));
  a->sigC   = malloc((cfg->M + 1) * sizeof(num_t));
  a->eps    = malloc((cfg->M + 1) * cfg->N * cfg->G * sizeof(num_t));
  a->eta    = malloc((cfg->M + 1) * cfg->G * sizeof(num_t));
  a->d      = malloc((cfg->M + 1) * sizeof(num_t));  
  a->tau    = malloc((cfg->M + 1) * sizeof(num_t));
  a->phi    = malloc((cfg->M + 1) * cfg->G * sizeof(num_t));
  a->thePhi = malloc((cfg->M + 1) * sizeof(num_t));
  a->sigPhi = malloc((cfg->M + 1) * sizeof(num_t));
  a->alp    = malloc((cfg->M + 1) * cfg->G * sizeof(num_t));
  a->theAlp = malloc((cfg->M + 1) * sizeof(num_t));
  a->sigAlp = malloc((cfg->M + 1) * sizeof(num_t));
  a->piAlp  = malloc((cfg->M + 1) * sizeof(num_t));
  a->del    = malloc((cfg->M + 1) * cfg->G * sizeof(num_t));
  a->theDel = malloc((cfg->M + 1) * sizeof(num_t));
  a->sigDel = malloc((cfg->M + 1) * sizeof(num_t));
  a->piDel  = malloc((cfg->M + 1) * sizeof(num_t));
  
  /* temporary and return values */
  
  a->tmp1 = malloc(cfg->G * sizeof(num_t));
  a->tmp2 = malloc(cfg->G * sizeof(num_t));

  a->Old   = malloc(cfg->N * sizeof(num_t));
  a->New   = malloc(cfg->N * sizeof(num_t));
  a->lOld  = malloc(cfg->N * sizeof(num_t));
  a->lNew  = malloc(cfg->N * sizeof(num_t));

  /* tuning parameters for Metropolis steps */
  
  a->tuneC = malloc(cfg->N * sizeof(num_t));
  a->tunePhi = malloc(cfg->G * sizeof(num_t));
  a->tuneEps = malloc(cfg->N * cfg->G * sizeof(num_t));

  /* number of acceptances for Metropolis steps */

  a->accC   = malloc(cfg->N * sizeof(int));
  a->accPhi = malloc(cfg->G * sizeof(int));
  a->accAlp = malloc(cfg->G * sizeof(int));
  a->accDel = malloc(cfg->G * sizeof(int));
  a->accEps = malloc(cfg->N * cfg->G * sizeof(int));
    
  return a;
}