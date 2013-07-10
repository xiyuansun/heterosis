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

__global__ void curand_setup_kernel(Chain *a, unsigned int seed){
  int G = a->G;
  int id = GENE;
  curand_init(seed, id, 0, &(a->states[id]));
}

__global__ void newChain_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int n, N = a->N, G = a->G;
  int g = GENE;
  num_t u;

  a->phi[iG(0, g)] = 1; /* rnormal(a->thePhi[0], a->sigPhi[0]);*/

  u = 1; /* runiform(0, 1) */;
  if(u < a->piAlp[0]){
    a->alp[iG(0, g)] = 0;
  } else {
    a->alp[iG(0, g)] = 1; /* rnormal(a->theAlp[0], a->sigAlp[0]); */
  }
    
  u = 1; /* runiform(0, 1); */
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

__host__ void newChain(Chain **host_a, Chain **dev_a, Config *cfg){ /* host */
  int n, g, G, *grp;
  count_t *y;
  num_t *lqts, s = 0, tmp, *tmpv, *yMeanG;

  y = readData(cfg);
  G = cfg->G;
  
  if(y == NULL)
    return;

  grp = readGrp(cfg);
  
  if(grp == NULL){
    free(y);
    return;
  }

  allocChainDevice(host_a, dev_a, cfg);
  
  /* data and configuration info */

  CUDA_CALL(cudaMemcpy(&((*dev_a)->M), &(cfg->M), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->N), &(cfg->N), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->G), &(cfg->G), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->burnin), &(cfg->burnin), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->heterosis), &(cfg->heterosis), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->someParmsFlag), &(cfg->someParmsFlag), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->allParmsFlag), &(cfg->allParmsFlag), sizeof(int), cudaMemcpyHostToDevice));  
  
  /* initialization constants */
  
  CUDA_CALL(cudaMemcpy(&((*dev_a)->sigC0), &(cfg->sigC0), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->d0), &(cfg->d0), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->aTau), &(cfg->aTau), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->aAlp), &(cfg->aAlp), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->aDel), &(cfg->aDel), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->bTau), &(cfg->bTau), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->bAlp), &(cfg->bAlp), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->bDel), &(cfg->bDel), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->gamPhi), &(cfg->gamPhi), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->gamAlp), &(cfg->gamAlp), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->gamDel), &(cfg->gamDel), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->sigPhi0), &(cfg->sigPhi0), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->sigAlp0), &(cfg->sigAlp0), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->sigDel0), &(cfg->sigDel0), sizeof(num_t), cudaMemcpyHostToDevice));
  
  /* choices to hold hyperparameters constant */
  
  CUDA_CALL(cudaMemcpy(&((*dev_a)->constSigC), &(cfg->constSigC), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->constD), &(cfg->constD), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->constTau), &(cfg->constTau), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->constThePhi), &(cfg->constThePhi), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->constTheAlp), &(cfg->constTheAlp), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->constTheDel), &(cfg->constTheDel), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->constSigPhi), &(cfg->constSigPhi), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->constSigAlp), &(cfg->constSigAlp), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->constSigDel), &(cfg->constSigDel), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->constPiAlp), &(cfg->constPiAlp), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->constPiDel), &(cfg->constPiDel), sizeof(int), cudaMemcpyHostToDevice));

  /* hyperparameters */
  
  CUDA_CALL(cudaMemcpy((*host_a)->sigC, &(cfg->sigC), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy((*host_a)->d, &(cfg->d), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy((*host_a)->tau, &(cfg->tau), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy((*host_a)->thePhi, &(cfg->thePhi), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy((*host_a)->theAlp, &(cfg->theAlp), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy((*host_a)->theDel, &(cfg->theDel), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy((*host_a)->sigPhi, &(cfg->sigPhi), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy((*host_a)->sigAlp, &(cfg->sigAlp), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy((*host_a)->sigDel, &(cfg->sigDel), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy((*host_a)->piAlp, &(cfg->piAlp), sizeof(num_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy((*host_a)->piDel, &(cfg->piDel), sizeof(num_t), cudaMemcpyHostToDevice));
  

  /* data */
  
  CUDA_CALL(cudaMemcpy((*host_a)->grp, grp, cfg->N * sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy((*host_a)->y, y, cfg->N * cfg->G * sizeof(int), cudaMemcpyHostToDevice));
  
  yMeanG = (num_t*) malloc(cfg->N * sizeof(num_t));
   
  for(n = 0; n < cfg->N; ++n){
    tmp = 0;
    
    for(g = 0; g < cfg->G; ++g)
      tmp += y[iG(n, g)];
    
    tmp /= cfg->G;
    yMeanG[n] = tmp;   
  }

  CUDA_CALL(cudaMemcpy((*host_a)->yMeanG, yMeanG, cfg->N * sizeof(num_t), cudaMemcpyHostToDevice));
  
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
  
  for(n = 0; n < cfg->N; ++n)
    tmpv[n] = lqts[n] - s;
  
  CUDA_CALL(cudaMemcpy((*host_a)->c, tmpv, cfg->N *sizeof(num_t), cudaMemcpyHostToDevice));
  
  /* set up curand states */
  
  curand_setup_kernel<<<NBLOCKS, NTHREADS>>>(*dev_a, cfg->seed);
  
  newChain_kernel1<<<NBLOCKS, NTHREADS>>>(*dev_a);
  newChain_kernel2<<<1, 1>>>(*dev_a);
  
  free(yMeanG);
  free(lqts);
  free(tmpv);
  free(grp);
  free(y);
}