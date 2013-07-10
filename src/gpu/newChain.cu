#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

__host__ int cmpfunc (const void *a, const void *b){
   return ( *(num_t*)a - *(num_t*)b );
}

__global__ void newChain_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int n, g = GENE, N = a->N, G = a->G;
  num_t u;

  a->phi[iG(0, g)] = 1; /* rnormal(a->thePhi[0], a->sigPhi[0]);*/

  u = 0.5; /* runiform(0, 1) */;
  if(u < a->piAlp[0]){
    a->alp[iG(0, g)] = 0;
  } else {
    a->alp[iG(0, g)] = 1; /* rnormal(a->theAlp[0], a->sigAlp[0]); */
  }
    
  u = 0.5; /* runiform(0, 1); */
  if(u < a->piDel[0]){
    a->del[iG(0, g)] = 0;
  } else {
    a->del[iG(0, g)] = 1; /* rnormal(a->theDel[0], a->sigDel[0]);*/
  }
 
  a->eta[iG(0, g)] = 1; /* 1/sqrt(rgamma(a->d[0] / 2, 
                 a->d[0] * a->tau[0] * a->tau[0] / 2, 0)); */

  for(n = 0; n < a->N; ++n)
    a->eps[iNG(0, n, g)] = 1; /* rnormal(0, a->eta[iG(0, g)]); */
    
}

__global__ void newChain_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  int n, g, G = a->G;

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

  a->tuneD = 100;

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

__host__ Chain *newChain(Config *cfg){ /* host */
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

  a = allocChain(cfg, 0);

  /* data and configuration info */

  CUDA_CALL(cudaMemcpy(&(a->M), &(cfg->M), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->N), &(cfg->N), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->G), &(cfg->G), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->burnin), &(cfg->burnin), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->heterosis), &(cfg->heterosis), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->someParmsFlag), &(cfg->someParmsFlag), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->allParmsFlag), &(cfg->allParmsFlag), sizeof(int), cudaMemcpyHostToDevice));  
  
  for(n = 0; n < cfg->N; ++n){
    printf("N = %d\n", cfg->N);
  
    CUDA_CALL(cudaMemcpy(&(a->grp[n]), &(grp[n]), sizeof(int), cudaMemcpyHostToDevice));
    tmp = 0;
    
    for(g = 0; g < cfg->G; ++g){
      a->y[iG(n, g)] = y[iG(n, g)];
      tmp += y[iG(n, g)];
    }
    tmp /= cfg->G;
    
    CUDA_CALL(cudaMemcpy(&(a->yMeanG[n]), &(tmp), sizeof(num_t), cudaMemcpyHostToDevice));
  }
    
    
    
  /* initialization constants */
  
  CUDA_CALL(cudaMemcpy(&(a->sigC0), &(cfg->sigC0), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->d0), &(cfg->d0), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->aTau), &(cfg->aTau), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->aAlp), &(cfg->aAlp), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->aDel), &(cfg->aDel), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->bTau), &(cfg->bTau), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->bAlp), &(cfg->bAlp), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->bDel), &(cfg->bDel), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->gamPhi), &(cfg->gamPhi), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->gamAlp), &(cfg->gamAlp), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->gamDel), &(cfg->gamDel), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->sigPhi0), &(cfg->sigPhi0), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->sigAlp0), &(cfg->sigAlp0), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->sigDel0), &(cfg->sigDel0), sizeof(num_t), cudaMemcpyHostToDevice));
  
  /* hyperparameters */
  
  CUDA_CALL(cudaMemcpy(&(a->sigC[0]), &(cfg->sigC), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->d[0]), &(cfg->d), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->tau[0]), &(cfg->tau), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->thePhi[0]), &(cfg->thePhi), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->theAlp[0]), &(cfg->theAlp), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->theDel[0]), &(cfg->theDel), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->sigPhi[0]), &(cfg->sigPhi), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->sigAlp[0]), &(cfg->sigAlp), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->sigDel[0]), &(cfg->sigDel), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->piAlp[0]), &(cfg->piAlp), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->piDel[0]), &(cfg->piDel), sizeof(num_t), cudaMemcpyHostToDevice));
  
  /* choices to hold hyperparameters constant */
  
  CUDA_CALL(cudaMemcpy(&(a->constSigC), &(cfg->constSigC), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->constD), &(cfg->constD), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->constTau), &(cfg->constTau), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->constThePhi), &(cfg->constThePhi), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->constTheAlp), &(cfg->constTheAlp), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->constTheDel), &(cfg->constTheDel), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->constSigPhi), &(cfg->constSigPhi), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->constSigAlp), &(cfg->constSigAlp), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->constSigDel), &(cfg->constSigDel), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->constPiAlp), &(cfg->constPiAlp), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&(a->constPiDel), &(cfg->constPiDel), sizeof(int), cudaMemcpyHostToDevice));
  
  lqts = (num_t*) malloc(cfg->N * sizeof(num_t));
  tmpv = (num_t*) malloc(cfg->G * sizeof(num_t));
  
  s = 0;
  for(n = 0; n < cfg->N; ++n){
    for(g = 0; g < cfg->G; ++g)
      tmpv[g] = y[iG(n, g)];
      
    qsort(tmpv, cfg->N, sizeof(num_t), cmpfunc);    
    lqts[n] = log(tmpv[(int) floor(cfg->G * 0.75)]);
    s += lqts[n];
  }
  
  s /= cfg->N;
  
  for(n = 0; n < cfg->N; ++n){
    tmp = lqts[n] - s;
    CUDA_CALL(cudaMemcpy(&(a->c[iN(0, n)]), &(tmp), sizeof(num_t), cudaMemcpyHostToDevice));
  }
  
  newChain_kernel1<<<NBLOCKS, NTHREADS>>>(a);
  newChain_kernel2<<<1, 1>>>(a);
  
  free(lqts);
  free(tmpv);
  free(grp);
  free(y);
    
  return a;
}