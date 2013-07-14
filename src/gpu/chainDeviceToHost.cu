#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

__host__ Chain *chainDeviceToHost(Chain *host_a, Chain *dev_a, Config *cfg){
  int N = cfg->N, G = cfg->G;
  Chain *allHost_a;
  
  if(cfg->verbose)
    printf("    Copying chain from device to host.\n");

  allHost_a = allocChainHost(cfg); 
  
  /* curand states */

  CUDA_CALL(cudaMemcpy(allHost_a->states, host_a->states, MAX_NG * sizeof(curandState_t), cudaMemcpyDeviceToHost));
  
  /* program options */

  CUDA_CALL(cudaMemcpy(&(allHost_a->m), &(dev_a->m), sizeof(int), cudaMemcpyDeviceToHost));  
  CUDA_CALL(cudaMemcpy(&(allHost_a->M), &(dev_a->M), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->N), &(dev_a->N), sizeof(int), cudaMemcpyDeviceToHost)); 
  CUDA_CALL(cudaMemcpy(&(allHost_a->G), &(dev_a->G), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->burnin), &(dev_a->burnin), sizeof(int), cudaMemcpyDeviceToHost));
  
  CUDA_CALL(cudaMemcpy(&(allHost_a->sumLogLik), &(dev_a->sumLogLik), sizeof(num_t), cudaMemcpyDeviceToHost));
  
  /* data */

  CUDA_CALL(cudaMemcpy(allHost_a->y, host_a->y, cfg->N * cfg->G * sizeof(count_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->yMeanG, host_a->yMeanG, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->grp, host_a->grp, cfg->N * sizeof(int), cudaMemcpyDeviceToHost));
  
  /* initialization constants */
  CUDA_CALL(cudaMemcpy(&(allHost_a->sigC0), &(dev_a->sigC0), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->d0), &(dev_a->d0), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->aTau), &(dev_a->aTau), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->aAlp), &(dev_a->aAlp), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->aDel), &(dev_a->aDel), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->bTau), &(dev_a->bTau), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->bAlp), &(dev_a->bAlp), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->bDel), &(dev_a->bDel), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->gamPhi), &(dev_a->gamPhi), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->gamAlp), &(dev_a->gamAlp), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->gamDel), &(dev_a->gamDel), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->sigPhi0), &(dev_a->sigPhi0), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->sigAlp0), &(dev_a->sigAlp0), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->sigDel0), &(dev_a->sigDel0), sizeof(num_t), cudaMemcpyDeviceToHost));
  
  /* parameters */
  
  CUDA_CALL(cudaMemcpy(allHost_a->c, host_a->c,  cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->sigC), &(dev_a->sigC),  sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->eps, host_a->eps,  cfg->N * cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->eta, host_a->eta,  cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->d), &(dev_a->d),  sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->tau), &(dev_a->tau),  sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->phi, host_a->phi,  cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->thePhi), &(dev_a->thePhi),  sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->sigPhi), &(dev_a->sigPhi),  sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->alp, host_a->alp,  cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->theAlp), &(dev_a->theAlp),  sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->sigAlp), &(dev_a->sigAlp),  sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->piAlp), &(dev_a->piAlp),  sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->del, host_a->del,  cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->theDel), &(dev_a->theDel),  sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->sigDel), &(dev_a->sigDel),  sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->piDel), &(dev_a->piDel),  sizeof(num_t), cudaMemcpyDeviceToHost));
  
  /* temporary and return values */
  
  CUDA_CALL(cudaMemcpy(&(allHost_a->s1), &(dev_a->s1), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->s2), &(dev_a->s2), sizeof(num_t), cudaMemcpyDeviceToHost));
  
  CUDA_CALL(cudaMemcpy(allHost_a->tmp1, host_a->tmp1, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->tmp2, host_a->tmp2, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  
  CUDA_CALL(cudaMemcpy(allHost_a->Old, host_a->Old, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->New, host_a->New, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->lOld, host_a->lOld, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->lNew, host_a->lNew, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
  
  /* current place in the chain */
    
  CUDA_CALL(cudaMemcpy(&(allHost_a->m), &(dev_a->m), sizeof(int), cudaMemcpyDeviceToHost));

  /* tuning parameters for Metropolis steps */
  
  CUDA_CALL(cudaMemcpy(&(allHost_a->tuneD), &(dev_a->tuneD), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->tuneC, host_a->tuneC, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->tunePhi, host_a->tunePhi, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->tuneEps, host_a->tuneEps, cfg->N * cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  
  /* number of acceptances for Metropolis steps */
  
  CUDA_CALL(cudaMemcpy(&(allHost_a->accD), &(dev_a->accD), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->accC, host_a->accC, cfg->N * sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->accPhi, host_a->accPhi, cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->accAlp, host_a->accAlp, cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->accDel, host_a->accDel, cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->accEps, host_a->accEps, cfg->N * cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
  
  /* choices to hold hyperparameters constant */ 
  
  CUDA_CALL(cudaMemcpy(&(allHost_a->constSigC), &(dev_a->constSigC), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->constD), &(dev_a->constD), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->constTau), &(dev_a->constTau), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->constThePhi), &(dev_a->constThePhi), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->constTheAlp), &(dev_a->constTheAlp), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->constTheDel), &(dev_a->constTheDel), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->constSigPhi), &(dev_a->constSigPhi), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->constSigAlp), &(dev_a->constSigAlp), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->constSigDel), &(dev_a->constSigDel), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->constPiAlp), &(dev_a->constPiAlp), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->constPiDel), &(dev_a->constPiDel), sizeof(int), cudaMemcpyDeviceToHost));

  CUDA_CALL(cudaMemcpy(allHost_a->dex, host_a->dex, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->hph, host_a->hph, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));  
  CUDA_CALL(cudaMemcpy(allHost_a->lph, host_a->lph, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->mph, host_a->mph, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  
  CUDA_CALL(cudaMemcpy(allHost_a->sumC, host_a->sumC, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->sumPhi, host_a->sumPhi, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));  
  CUDA_CALL(cudaMemcpy(allHost_a->sumAlp, host_a->sumAlp, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->sumDel, host_a->sumDel, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->sumEps, host_a->sumEps, cfg->N * cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  
  return allHost_a;
}