#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <cuda.h>
#include <curand_kernel.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

__host__ int cmpfunc (const void *a, const void *b){
   return ( *(num_t*)a - *(num_t*)b );
}

__global__ void curand_setup_kernel(Chain *a, int *seeds){ /* kernel <<<G, 1>>> */
  int id = IDX, N = a->N, G = a->G;
  if(id < MAX_NG)
    curand_init(seeds[id], id, 0, &(a->states[id]));
}

__host__ void copyHyperParms(Chain *host_a, Chain *dev_a, Config *cfg){

  if(cfg->constSigC)
    CUDA_CALL(cudaMemcpy(&((*dev_a)->sigC), &(cfg->sigC), sizeof(num_t), cudaMemcpyHostToDevice));

  if(cfg->constD)
    CUDA_CALL(cudaMemcpy(&((*dev_a)->d), &(cfg->d), sizeof(num_t), cudaMemcpyHostToDevice));

  if(cfg->constTau)
    CUDA_CALL(cudaMemcpy(&((*dev_a)->tau), &(cfg->tau), sizeof(num_t), cudaMemcpyHostToDevice));

  if(cfg->constThePhi)
    CUDA_CALL(cudaMemcpy(&((*dev_a)->thePhi), &(cfg->thePhi), sizeof(num_t), cudaMemcpyHostToDevice));

  if(cfg->constTheAlp)
    CUDA_CALL(cudaMemcpy(&((*dev_a)->theAlp), &(cfg->theAlp), sizeof(num_t), cudaMemcpyHostToDevice));

  if(cfg->constTheDel)
    CUDA_CALL(cudaMemcpy(&((*dev_a)->theDel), &(cfg->theDel), sizeof(num_t), cudaMemcpyHostToDevice));

  if(cfg->constSigPhi)
    CUDA_CALL(cudaMemcpy(&((*dev_a)->sigPhi), &(cfg->sigPhi), sizeof(num_t), cudaMemcpyHostToDevice));

  if(cfg->constSigAlp)
    CUDA_CALL(cudaMemcpy(&((*dev_a)->sigAlp), &(cfg->sigAlp), sizeof(num_t), cudaMemcpyHostToDevice));

  if(cfg->constSigDel)
    CUDA_CALL(cudaMemcpy(&((*dev_a)->sigDel), &(cfg->sigDel), sizeof(num_t), cudaMemcpyHostToDevice));

  if(cfg->constPiAlp)  
    CUDA_CALL(cudaMemcpy(&((*dev_a)->piAlp), &(cfg->piAlp), sizeof(num_t), cudaMemcpyHostToDevice));

  if(cfg->constPiDel)
    CUDA_CALL(cudaMemcpy(&((*dev_a)->piDel), &(cfg->piDel), sizeof(num_t), cudaMemcpyHostToDevice));
}

__global__ void newChain_kernel1(Chain *a){ /* kernel <<<1, 1>>> */
  int n;

  a->m = 1; 
  a->accD = 0;
  a->tuneD = 400;
  
  a->meanLogLik = 0;
  a->logLikMean = 0;
  a->dic = 0;
  
  for(n = 0; n < a->N; ++n){
    a->c[n] = 0;
    a->accC[n] = 0;
    a->tuneC[n] = 1;
  }
  
  if(!a->constTau)
    a->tau = sqrt(rgammaDevice(a, 1, a->aTau, a->bTau, 0));
 
  if(!a->constPiAlp)
    a->piAlp = rbetaDevice(a, 1, a->aAlp, a->bAlp);
  
  if(!a->constPiDel)
    a->piDel = rbetaDevice(a, 1, a->aDel, a->bDel);

  if(!a->constD)
    a->d = runiformDevice(a, 1, 0, a->d0);
 
  if(!a->constThePhi)
    a->thePhi = rnormalDevice(a, 1, 0, a->gamPhi);

  if(!a->constTheAlp)
    a->theAlp = rnormalDevice(a, 1, 0, a->gamAlp);

  if(!a->constTheDel)
    a->theDel = rnormalDevice(a, 1, 0, a->gamDel);
 
  if(!a->constSigC)
    a->sigC = runiformDevice(a, 1, 0, a->sigC0);
   
  if(!a->constSigPhi)
    a->sigPhi = runiformDevice(a, 1, 0, a->sigPhi0);

  if(!a->constSigAlp)
    a->sigAlp = runiformDevice(a, 1, 0, a->sigAlp0);
 
  if(!a->constSigDel)
    a->sigDel = runiformDevice(a, 1, 0, a->sigDel0);  
}

__global__ void newChain_kernel2(Chain *a){ /* kernel <<<G, 1>>> */
  int n, g = IDX, G = a->G;
  num_t u;

  if(g < G){
    a->dex[g] = 0;
    a->hph[g] = 0;
    a->lph[g] = 0;
    a->mph[g] = 0;

    a->phi[g] = rnormalDevice(a, g, a->thePhi, a->sigPhi);
    a->eta[g] = 1/sqrt(rgammaDevice(a, g, a->d / 2, 
                   a->d * a->tau * a->tau / 2, 0));

    a->accPhi[g] = 0;
    a->accAlp[g] = 0;
    a->accDel[g] = 0;

    a->tunePhi[g] = 1; 

    for(n = 0; n < a->N; ++n){
      a->accEps[iG(n, g)] = 0;
      a->tuneEps[iG(n, g)] = 1;
      a->eps[iG(n, g)] = rnormalDevice(a, g, 0, a->eta[g]);
    }
    
    u = runiformDevice(a, g, 0, 1);
    if(u < a->piAlp){
      a->alp[g] = 0;
    } else {
      a->alp[g] = rnormalDevice(a, g, a->theAlp, a->sigAlp);
    }
    
    u = runiformDevice(a, g, 0, 1);
    if(u < a->piDel){
      a->del[g] = 0;
    } else {
      a->del[g] = rnormalDevice(a, g, a->theDel, a->sigDel);
    }
  }
}

void newChain(Chain **host_a, Chain **dev_a, Config *cfg){ /* host */
  int n, g, N, G, i, *grp, *seeds, *dev_seeds;
  count_t *y;
  num_t tmp, *yMeanG;

  y = readData(cfg);
  
  N = cfg->N;
  G = cfg->G;
  
  if(y == NULL)
    return;

  grp = readGrp(cfg);
  
  if(grp == NULL){
    free(y);
    return;
  }
  
  if(cfg->verbose)
    printf("  Allocating chain.\n"); 
    
  allocChainDevice(host_a, dev_a, cfg);
  
  /* configuration info */
  
  CUDA_CALL(cudaMemcpy(&((*dev_a)->M), &(cfg->M), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->N), &(cfg->N), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->G), &(cfg->G), sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(&((*dev_a)->burnin), &(cfg->burnin), sizeof(int), cudaMemcpyHostToDevice)); 
  CUDA_CALL(cudaMemcpy(&((*dev_a)->heterosis), &(cfg->heterosis), sizeof(int), cudaMemcpyHostToDevice)); 
    
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
  
  copyHyperParms(host_a, dev_a, cfg);
    
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
  
  /* set up CURAND */ 
  
  seeds = (int*) malloc(MAX_NG * sizeof(int));
  CUDA_CALL(cudaMalloc((void**) &dev_seeds, MAX_NG * sizeof(int)));  

  for(i = 0; i < MAX_NG; ++i)
    seeds[i] = rand();
    
  CUDA_CALL(cudaMemcpy(dev_seeds, seeds, MAX_NG * sizeof(int), cudaMemcpyHostToDevice));
  curand_setup_kernel<<<NG_GRID, NG_BLOCK>>>(*dev_a, dev_seeds);
  
  /* set up the bulk of the chain */
  
  newChain_kernel1<<<1, 1>>>(*dev_a);
  newChain_kernel2<<<G_GRID, G_BLOCK>>>(*dev_a);
 
  free(yMeanG);
  free(grp);
  free(y);
  free(seeds);
  cudaFree(dev_seeds);
}