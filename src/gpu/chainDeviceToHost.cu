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

  Chain *allHost_a;
  
  float myTime;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  
  fprintf(cfg->log, "  Copying chain from device to host. Takes a LONG time.\n");

  allocChainHost(&allHost_a, cfg); 
  
  /* program options */
  
  CUDA_CALL(cudaMemcpy(&(allHost_a->M), &(dev_a->M), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->N), &(dev_a->N), sizeof(int), cudaMemcpyDeviceToHost)); 
  CUDA_CALL(cudaMemcpy(&(allHost_a->G), &(dev_a->G), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->burnin), &(dev_a->burnin), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->heterosis), &(dev_a->heterosis), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->parmsFlag), &(dev_a->parmsFlag), sizeof(int), cudaMemcpyDeviceToHost));
  
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
  
  CUDA_CALL(cudaMemcpy(allHost_a->c, host_a->c, (cfg->M + 1) * cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->sigC, host_a->sigC, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->eps, host_a->eps, (cfg->M + 1) * cfg->N * cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->eta, host_a->eta, (cfg->M + 1) * cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->d, host_a->d, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->tau, host_a->tau, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->phi, host_a->phi, (cfg->M + 1) * cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->thePhi, host_a->thePhi, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->sigPhi, host_a->sigPhi, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->alp, host_a->alp, (cfg->M + 1) * cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->theAlp, host_a->theAlp, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->sigAlp, host_a->sigAlp, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->piAlp, host_a->piAlp, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->del, host_a->del, (cfg->M + 1) * cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->theDel, host_a->theDel, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->sigDel, host_a->sigDel, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->piDel, host_a->piDel, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  
  /* temporary and return values */
  
  CUDA_CALL(cudaMemcpy(&(allHost_a->s1), &(dev_a->s1), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->s2), &(dev_a->s2), sizeof(num_t), cudaMemcpyDeviceToHost));
  
  CUDA_CALL(cudaMemcpy(allHost_a->tmp1, host_a->tmp1, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->tmp2, host_a->tmp2, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  
  CUDA_CALL(cudaMemcpy(allHost_a->Old, host_a->Old, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->New, host_a->New, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->lOld, host_a->lOld, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(allHost_a->lNew, host_a->lNew, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
  
  /* current place in the chain of each parameter */
    
  CUDA_CALL(cudaMemcpy(&(allHost_a->mC), &(dev_a->mC), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->mSigC), &(dev_a->mSigC), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->mEps), &(dev_a->mEps), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->mEta), &(dev_a->mEta), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->mD), &(dev_a->mD), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->mTau), &(dev_a->mTau), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->mPhi), &(dev_a->mPhi), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->mAlp), &(dev_a->mAlp), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->mDel), &(dev_a->mDel), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->mThePhi), &(dev_a->mThePhi), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->mTheAlp), &(dev_a->mTheAlp), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->mTheDel), &(dev_a->mTheDel), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->mSigPhi), &(dev_a->mSigPhi), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->mSigAlp), &(dev_a->mSigAlp), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->mSigDel), &(dev_a->mSigDel), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->mPiAlp), &(dev_a->mPiAlp), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(allHost_a->mPiDel), &(dev_a->mPiDel), sizeof(int), cudaMemcpyDeviceToHost));
  
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

  /* curand states */

  CUDA_CALL(cudaMemcpy(allHost_a->states, host_a->states, cfg->G * sizeof(curandState), cudaMemcpyDeviceToHost));

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  
  printf("copy\n");
  
  fprintf(cfg->time, "%0.3f ", myTime); /* elapsed time in minutes */

  return allHost_a;
}