#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int cmpfunc (const void *a, const void *b){
   return ( *(num_t*)a - *(num_t*)b );
}

void newChain_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int n, g, G = a->G;
  num_t u;

  for(g = 0; g < a->G; ++g){
    a->dex[g] = 0;
    a->hph[g] = 0;
    a->lph[g] = 0;
    a->mph[g] = 0;

    a->phi[g] = rnormal(a->thePhi, a->sigPhi);

    u = runiform(0, 1);
    if(u < a->piAlp){
      a->alp[g] = 0;
    } else {
      a->alp[g] = rnormal(a->theAlp, a->sigAlp);
    }
    
    u = runiform(0, 1);
    if(u < a->piDel){
      a->del[g] = 0;
    } else {
      a->del[g] = rnormal(a->theDel, a->sigDel);
    }
 
    a->eta[g] = 1/sqrt(rgamma(a->d / 2, 
                   a->d * a->tau * a->tau / 2, 0));

    for(n = 0; n < a->N; ++n)
      a->eps[iG(n, g)] = rnormal(0, a->eta[g]);
    
  }
}

void newChain_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  int n, g, G = a->G;

  a->m = 1;
  
  /* counts toward differential expression and heterosis */
  
  a->tuneD = 400;
  
  for(n = 0; n < a->N; ++n)
    a->tuneC[n] = 1;

  for(g = 0; g < a->G; ++g){
    a->tunePhi[g] = 1;

    for(n = 0; n < a->N; ++n)
      a->tuneEps[iG(n, g)] = 1;
  }
  
  a->accD = 0;

  for(n = 0; n < a->N; ++n){
    a->accC[n] = 0;
  
    for(g = 0; g < a->G; ++g)
      a->accEps[iG(n, g)] = 0;
  }

  for(g = 0; g < a->G; ++g){
    a->accPhi[g] = 0;
    a->accAlp[g] = 0;
    a->accDel[g] = 0;
  }
}

Chain *newChain(Config *cfg){ /* host */
  int n, g, N, G, *grp;
  count_t *y;
  num_t *lqts, s = 0, tmp, *tmpv;
  Chain *a;

  y = readData(cfg);
  
  N = cfg->N;
  G = cfg->G;
  
  if(y == NULL)
    return NULL;

  grp = readGrp(cfg);
  
  if(grp == NULL){
    free(y);
    
    return NULL;
  }
  
  a = allocChain(cfg);

  /* data and configuration info */
  
  a->M = cfg->M;
  a->N = cfg->N;
  a->G = cfg->G;
  a->burnin = cfg->burnin;
  a->heterosis = cfg->heterosis;
  a->parmsFlag = cfg->parmsFlag;
  
  for(n = 0; n < cfg->N; ++n){
    
    a->grp[n] = grp[n];
    tmp = 0;
    
    for(g = 0; g < cfg->G; ++g){
      a->y[iG(n, g)] = y[iG(n, g)];
      tmp += y[iG(n, g)];
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
  
  a->sigC   = cfg->sigC;
  a->d      = cfg->d;
  a->tau    = cfg->tau;
  a->thePhi = cfg->thePhi;
  a->theAlp = cfg->theAlp;
  a->theDel = cfg->theDel;
  a->sigPhi = cfg->sigPhi;
  a->sigAlp = cfg->sigAlp;
  a->sigDel = cfg->sigDel;
  a->piAlp  = cfg->piAlp;
  a->piDel  = cfg->piDel;
  
  /* choices to hold hyperparameters constant */
  
  a->constSigC   = cfg->constSigC;
  a->constD      = cfg->constD;
  a->constTau    = cfg->constTau;
  a->constThePhi = cfg->constThePhi;
  a->constTheAlp = cfg->constTheAlp;
  a->constTheDel = cfg->constTheDel;
  a->constSigPhi = cfg->constSigPhi;
  a->constSigAlp = cfg->constSigAlp;
  a->constSigDel = cfg->constSigDel;
  a->constPiAlp  = cfg->constPiAlp;
  a->constPiDel  = cfg->constPiDel;
  
  /* initial normalization factors, c */
  
  lqts = (num_t*) malloc(cfg->N * sizeof(num_t));
  tmpv = (num_t*) malloc(cfg->G * sizeof(num_t));
  
  s = 0;
  for(n = 0; n < cfg->N; ++n){
    for(g = 0; g < cfg->G; ++g)
      tmpv[g] = y[iG(n, g)];
      
    qsort(tmpv, cfg->G, sizeof(num_t), cmpfunc);   
     
    lqts[n] = log(tmpv[(int) floor(cfg->G * 0.75)]);
    s += lqts[n];
  }
  
  s /= cfg->N;
  
  for(n = 0; n < cfg->N; ++n)
    a->c[n] = lqts[n] - s;
  
  newChain_kernel1(a);
  newChain_kernel2(a);
  
  free(lqts);
  free(tmpv);
  free(grp);
  free(y);

  return a;
}