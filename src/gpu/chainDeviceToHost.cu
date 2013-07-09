#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

__host__ Chain *chainDeviceToHost(Chain *a, Config *cfg){

  Chain *host_a = allocChain(cfg, 1); 
  
  /* program options */
  
  CUDA_CALL(cudaMemcpy(&(host_a->M), &(a->M), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->N), &(a->N), sizeof(int), cudaMemcpyDeviceToHost)); 
  CUDA_CALL(cudaMemcpy(&(host_a->G), &(a->G), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->burnin), &(a->burnin), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->heterosis), &(a->heterosis), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->someParmsFlag), &(a->someParmsFlag), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->allParmsFlag), &(a->allParmsFlag), sizeof(int), cudaMemcpyDeviceToHost));
  
  /* data */
printf("2\n");
  CUDA_CALL(cudaMemcpy(host_a->y, a->y, cfg->N * cfg->G * sizeof(count_t), cudaMemcpyDeviceToHost)); printf("3\n"); 
  CUDA_CALL(cudaMemcpy(host_a->yMeanG, a->yMeanG, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));printf("4\n");
  CUDA_CALL(cudaMemcpy(host_a->grp, a->grp, cfg->N * sizeof(int), cudaMemcpyDeviceToHost));printf("5\n");
  
  /* initialization constants */
  CUDA_CALL(cudaMemcpy(&(host_a->sigC0), &(a->sigC0), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->d0), &(a->d0), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->aTau), &(a->aTau), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->aAlp), &(a->aAlp), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->aDel), &(a->aDel), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->bTau), &(a->bTau), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->bAlp), &(a->bAlp), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->bDel), &(a->bDel), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->gamPhi), &(a->gamPhi), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->gamAlp), &(a->gamAlp), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->gamDel), &(a->gamDel), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->sigPhi0), &(a->sigPhi0), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->sigAlp0), &(a->sigAlp0), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->sigDel0), &(a->sigDel0), sizeof(num_t), cudaMemcpyDeviceToHost));
  
  /* parameters */
  
  CUDA_CALL(cudaMemcpy(host_a->c, a->c, (cfg->M + 1) * cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->sigC, a->sigC, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->eps, a->eps, (cfg->M + 1) * cfg->N * cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->eta, a->eta, (cfg->M + 1) * cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->d, a->d, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->tau, a->tau, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->phi, a->phi, (cfg->M + 1) * cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->thePhi, a->thePhi, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->sigPhi, a->sigPhi, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->alp, a->alp, (cfg->M + 1) * cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->theAlp, a->theAlp, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->sigAlp, a->sigAlp, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->piAlp, a->piAlp, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->del, a->del, (cfg->M + 1) * cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->theDel, a->theDel, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->sigDel, a->sigDel, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->piDel, a->piDel, (cfg->M + 1) * sizeof(num_t), cudaMemcpyDeviceToHost));
  
  /* temporary and return values */
  
  CUDA_CALL(cudaMemcpy(&(host_a->s1), &(a->s1), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->s2), &(a->s2), sizeof(num_t), cudaMemcpyDeviceToHost));
  
  CUDA_CALL(cudaMemcpy(host_a->tmp1, a->tmp1, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->tmp2, a->tmp2, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  
  CUDA_CALL(cudaMemcpy(host_a->Old, a->Old, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->New, a->New, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->lOld, a->lOld, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->lNew, a->lNew, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
  
  /* current place in the chain of each parameter */
    
  CUDA_CALL(cudaMemcpy(&(host_a->mC), &(a->mC), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->mSigC), &(a->mSigC), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->mEps), &(a->mEps), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->mEta), &(a->mEta), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->mD), &(a->mD), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->mTau), &(a->mTau), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->mPhi), &(a->mPhi), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->mAlp), &(a->mAlp), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->mDel), &(a->mDel), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->mThePhi), &(a->mThePhi), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->mTheAlp), &(a->mTheAlp), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->mTheDel), &(a->mTheDel), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->mSigPhi), &(a->mSigPhi), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->mSigAlp), &(a->mSigAlp), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->mSigDel), &(a->mSigDel), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->mPiAlp), &(a->mPiAlp), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->mPiDel), &(a->mPiDel), sizeof(int), cudaMemcpyDeviceToHost));
  
  /* tuning parameters for Metropolis steps */
  
  CUDA_CALL(cudaMemcpy(&(host_a->tuneD), &(a->tuneD), sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->tuneC, a->tuneC, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->tunePhi, a->tunePhi, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->tuneEps, a->tuneEps, cfg->N * cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  
  /* number of acceptances for Metropolis steps */
  
  CUDA_CALL(cudaMemcpy(&(host_a->accD), &(a->accD), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->accC, a->accC, cfg->N * sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->accPhi, a->accPhi, cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->accAlp, a->accAlp, cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->accDel, a->accDel, cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(host_a->accEps, a->accEps, cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
  
  /* choices to hold hyperparameters constant */
  
  CUDA_CALL(cudaMemcpy(&(host_a->constSigC), &(a->constSigC), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->constD), &(a->constD), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->constTau), &(a->constTau), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->constThePhi), &(a->constThePhi), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->constTheAlp), &(a->constTheAlp), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->constTheDel), &(a->constTheDel), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->constSigPhi), &(a->constSigPhi), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->constSigAlp), &(a->constSigAlp), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->constSigDel), &(a->constSigDel), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->constPiAlp), &(a->constPiAlp), sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(&(host_a->constPiDel), &(a->constPiDel), sizeof(int), cudaMemcpyDeviceToHost));
  
  return a;
}