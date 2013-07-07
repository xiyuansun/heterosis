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

void newChainHost_kernel1(Chain *a){
  int n, g;
  num_t u;

  for(g = 0; g < a->G; ++g){
    a->phi[0][g] = normalHost(a->thePhi[0], a->sigPhi[0]);

    u = uniformHost(0, 1);
    if(u < a->piAlp[0]){
      a->alp[0][g] = 0;
    } else {
      a->alp[0][g] = normalHost(a->theAlp[0], a->sigAlp[0]);
    }
    
    u = uniformHost(0, 1);
    if(u < a->piDel[0]){
      a->del[0][g] = 0;
    } else {
      a->del[0][g] = normalHost(a->theDel[0], a->sigDel[0]);
    }
 
    a->eta[0][g] = 1/sqrt(gammaHost(a->d[0] / 2, 
                   a->d[0] * a->tau[0] * a->tau[0] / 2, 0));

    for(n = 0; n < a->N; ++n)
      a->eps[0][n][g] = normalHost(0, a->eta[0][g]);
  }
}

void newChainHost_kernel2(Chain *a){
  int n, g;

  a->mC = 0;
  a->mSigC = 0;

  a->mEps = 0;
  a->mEta = 0;
  a->mD = 0;
  a->mTau = 0;

  a->mPhi = 0;
  a->mAlp = 0;
  a->mDel = 0;

  a->mThePhi = 0;
  a->mTheAlp = 0;
  a->mTheDel = 0;

  a->mSigPhi = 0;
  a->mSigAlp = 0;
  a->mSigDel = 0;

  a->mPiAlp = 0;
  a->mPiDel = 0;

  a->tuneD = 500;

  for(n = 0; n < a->N; ++n)
    a->tuneC[n] = 1;

  for(g = 0; g < a->G; ++g){
    a->tunePhi[g] = 1;

    for(n = 0; n < a->N; ++n)
      a->tuneEps[n][g] = 1;
  }
  
  a->accD = 0;

  for(n = 0; n < a->N; ++n){
    a->accC[n] = 0;
  
    for(g = 0; g < a->G; ++g)
      a->accEps[n][g] = 0;
  }

  for(g = 0; g < a->G; ++g){
    a->accPhi[g] = 0;
    a->accAlp[g] = 0;
    a->accDel[g] = 0;
  }
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
  
  newChainHost_kernel1(a);
  newChainHost_kernel2(a);
    
  return a;
}