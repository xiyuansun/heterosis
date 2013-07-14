#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <curand_kernel.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

__global__ void curand_setup_kernel(curandState_t *states, int *seeds, int N, int G){ /* kernel <<<G, 1>>> */
  int id = IDX;
  if(id < MAX_NG)
    curand_init(seeds[id], id, 0, &(states[id]));
}

Config *config(int argc, char **argv){
  int n, g, i, N, G, *seeds, *dev_seeds;
  num_t tmp;
  curandState_t *states;
  
  Config *cfg = (Config*) malloc(sizeof(Config));
  cfg->chainNum = 1;
  
  /* default filenames */        

  strcpy(cfg->dataFile, "../data/data.txt"); 
  strcpy(cfg->groupFile, "../data/grp.txt");
   
  cfg->ratesFlag = 0;
  cfg->hyperFlag = 0;
  cfg->parmsFlag = 0;
  cfg->timeFlag = 0;
  cfg->verbose = 0;
  cfg->diagnostics = 0;

  cfg->chains = 2;
  cfg->m = 1;
  cfg->M = 10;
  cfg->N = 0;
  cfg->G = 0;
  cfg->burnin = cfg->M / 2;
  cfg->joint = 0;
  cfg->seed = 22;

  /* default initialization constants */

  cfg->sigC0 = 10;
  cfg->d0 = 1000;
  cfg->aTau = 100;
  cfg->aAlp = 1;
  cfg->aDel = 1;
  cfg->bTau = 100;
  cfg->bAlp = 1;
  cfg->bDel = 1;
  cfg->gamPhi = 2;
  cfg->gamAlp = 2;
  cfg->gamDel = 2;
  cfg->sigPhi0 = 2;
  cfg->sigAlp0 = 2;
  cfg->sigDel0 = 2;
  
  /* default: don't hold hyperparameters constant */
  
  cfg->constSigC = 0;
  cfg->constD = 0;
  cfg->constTau = 0;
  cfg->constThePhi = 0;
  cfg->constTheAlp = 0;
  cfg->constTheDel = 0;
  cfg->constSigPhi = 0;
  cfg->constSigAlp = 0;
  cfg->constSigDel = 0;
  cfg->constPiAlp = 0;
  cfg->constPiDel = 0;
  
  cfg->timeC = 0;
  cfg->timeTau = 0;
  cfg->timePiAlp = 0;
  cfg->timePiDel = 0;
  cfg->timeD = 0;
  cfg->timeThePhi = 0;
  cfg->timeTheAlp = 0;
  cfg->timeTheDel = 0;
  cfg->timeSigC = 0;
  cfg->timeSigPhi = 0;
  cfg->timeSigAlp = 0;
  cfg->timeSigDel = 0;
  cfg->timeEta = 0;
  cfg->timeEps = 0;
  cfg->timePhi = 0;
  cfg->timeAlp = 0;
  cfg->timeDel = 0;

  getopts(cfg, argc, argv);
  srand(cfg->seed);
  
  /* read data and grp */
  
  cfg->y = readData(cfg);
  
  N = cfg->N;
  G = cfg->G;
  
  if(cfg->y == NULL){
    return NULL;
  }

printf("outside readdata %d %d\n", cfg->N, cfg->G);
  cfg->grp = readGrp(cfg);
  
  
  if(cfg->grp == NULL){
    free(cfg->y);
    return NULL;
  }
  
  allocConfig(cfg);

  for(n = 0; n < cfg->N; ++n){
    tmp = 0;
    
    for(g = 0; g < cfg->G; ++g)
      tmp += cfg->y[iG(n, g)];
    
    tmp /= cfg->G;
    cfg->yMeanG[n] = tmp;   
  }
  
  CUDA_CALL(cudaMemcpy(cfg->devY, cfg->y, cfg->N * cfg->G * sizeof(count_t), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(cfg->devGrp, cfg->grp, cfg->N * sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(cfg->devYMeanG, cfg->yMeanG, cfg->N * sizeof(int), cudaMemcpyHostToDevice));  
    
  /* set up CURAND */
  
  seeds = (int*) malloc(MAX_NG * sizeof(int));
  CUDA_CALL(cudaMalloc((void**) &dev_seeds, MAX_NG * sizeof(int)));  
  CUDA_CALL(cudaMalloc((void**) &states, MAX_NG * sizeof(curandState_t)));  

  for(i = 0; i < MAX_NG; ++i)
    seeds[i] = rand(); 
    
  CUDA_CALL(cudaMemcpy(dev_seeds, seeds, MAX_NG * sizeof(int), cudaMemcpyHostToDevice));
  curand_setup_kernel<<<NG_GRID, NG_BLOCK>>>(states, dev_seeds, cfg->N, cfg->G);
  cfg->states = states;

  /* 
   *  All hyperparameters set in getopts() will be treated as constant.
   *  All the others must be given initial values.
   */
   
  if(!cfg->constTau)
    cfg->tau = sqrt(rgamma(cfg->aTau, cfg->bTau, 0));
 
  if(!cfg->constPiAlp)
    cfg->piAlp = rbeta(cfg->aAlp, cfg->bAlp);
  
  if(!cfg->constPiDel)
    cfg->piDel = rbeta(cfg->aDel, cfg->bDel);

  if(!cfg->constD)
    cfg->d = runiform(0, cfg->d0);
 
  if(!cfg->constThePhi)
    cfg->thePhi = rnormal(0, cfg->gamPhi);

  if(!cfg->constTheAlp)
    cfg->theAlp = rnormal(0, cfg->gamAlp);

  if(!cfg->constTheDel)
    cfg->theDel = rnormal(0, cfg->gamDel);
 
  if(!cfg->constSigC)
    cfg->sigC = runiform(0, cfg->sigC0);
  
  if(!cfg->constSigPhi)
    cfg->sigPhi = runiform(0, cfg->sigPhi0);

  if(!cfg->constSigAlp)
    cfg->sigAlp = runiform(0, cfg->sigAlp0);

  if(!cfg->constSigDel)
    cfg->sigDel = runiform(0, cfg->sigDel0);
   
  system("mkdir -p ../out/");
  system("mkdir -p ../out/probs/");
  
  if(cfg->ratesFlag)
    system("mkdir -p ../out/rates/");
  
  if(cfg->hyperFlag)
    system("mkdir -p ../out/hyper/");
  
  if(cfg->parmsFlag)
    system("mkdir -p ../out/parms/"); 
  
  if(cfg->timeFlag)
    system("mkdir -p ../out/time/");
  
  if(cfg->diagnostics){
    system("mkdir -p ../out/diagnostics/");
    system("rm -f ../out/diagnostics/dic.txt");
  }
  
  free(seeds);
  cudaFree(dev_seeds);
  return cfg;
}