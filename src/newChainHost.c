#include <Chain.h>
#include <Config.h>
#include <functions.h>
#include <math.h>
#include <numericTypes.h>
#include <stdio.h>
#include <stdlib.h>

int cmpfunc (const void *a, const void *b){
   return ( *(num_t*)a - *(num_t*)b );
}

Chain *newChainHost(Config *cfg){
  int n, g, *grp;
  count_t **y;
  num_t *lqts, s = 0, tmp, *tmpv;
  Chain *a;

  y = readData(cfg);
  
  if(y == NULL)
    return NULL;

  grp = readGrp(cfg);
  
  if(grp == NULL){
    for(n = 0; n < cfg->N; ++n)
      free(y[n]);
    free(y);
    
    return NULL;
  }

  a = allocChainHost(cfg);

  /* data and configuration info */
  
  a->M = cfg->M;
  a->N = cfg->N;
  a->G = cfg->G;
  a->burnin = cfg->burnin;
  a->heterosis = cfg->heterosis;
  a->someParmsFlag = cfg->someParmsFlag;
  a->allParmsFlag = cfg->allParmsFlag;
  
  for(n = 0; n < cfg->N; ++n){
    a->grp[n] = grp[n];
    tmp = 0;
    
    for(g = 0; g < cfg->G; ++g){
      a->y[n][g] = y[n][g];
      tmp += y[n][g];
    }
           
    a->yMeanG[n] = tmp / cfg->G;
  }
    
  /* initialization constants */
  
  a->sigC0   = cfg->sigC0;
  a->d0      = cfg->d0;
  a->aTau    = cfg->aTau;
  a->aAlp    = cfg->aAlp;
  a->aDel    = cfg->aDel; 
  a->bTau    = cfg->bTau;
  a->bAlp    = cfg->bAlp;
  a->bDel    = cfg->bDel;  
  a->gamPhi  = cfg->gamPhi;
  a->gamAlp  = cfg->gamAlp;
  a->gamDel  = cfg->gamDel;
  a->sigPhi0 = cfg->sigPhi0;
  a->sigAlp0 = cfg->sigAlp0;
  a->sigDel0 = cfg->sigDel0;
  
  /* hyperparameters */
  
  a->sigC[0]   = cfg->sigC;
  a->d[0]      = cfg->d;
  a->tau[0]    = cfg->tau;
  a->thePhi[0] = cfg->thePhi;
  a->theAlp[0] = cfg->theAlp;
  a->theDel[0] = cfg->theDel;
  a->sigPhi[0] = cfg->sigPhi;
  a->sigAlp[0] = cfg->sigAlp;
  a->sigDel[0] = cfg->sigDel;
  a->piAlp[0]  = cfg->piAlp;
  a->piDel[0]  = cfg->piDel;
  
  lqts = malloc(cfg->N * sizeof(num_t));
  tmpv = malloc(cfg->G * sizeof(num_t));
  
  s = 0;
  for(n = 0; n < cfg->N; ++n){
    for(g = 0; g < cfg->G; ++g)
      tmpv[g] = y[n][g];
      
    qsort(tmpv, cfg->N, sizeof(num_t), cmpfunc);    
    lqts[n] = log(tmpv[(int) floor(cfg->G * 0.75)]);
    s += lqts[n];
  }
  
  s /= cfg->N;
  
  for(n = 0; n < cfg->N; ++n)
    a->c[0][n] = lqts[n] - s;
  
  return a;
}