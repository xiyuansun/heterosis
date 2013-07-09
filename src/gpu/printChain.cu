#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

__host__ void printChain(Chain *a, Config *cfg, int onHost){

  Chain *host_a = onHost ? a : chainDeviceToHost(a, cfg);

  printf("M = %d\n", host_a->M); ;
  printf("N = %d\n", host_a->N);
  printf("G = %d\n\n", host_a->G);  
  
  printf("burnin = %d\n", host_a->burnin);
  printf("heterosis = %d\n", host_a->heterosis);
  printf("someParmsFlag = %d\n", host_a->someParmsFlag);
  printf("allParmsFlag = %d\n\n", host_a->allParmsFlag); 

  pi2(host_a->y, host_a->N, host_a->G, "y = \n");
  pf1(host_a->yMeanG, host_a->N, "yMeanG =\n");
  pi1(host_a->grp, host_a->N, "grp =\n");
  
  printf("sigC0   = "); printf(NUM_TF, host_a->sigC0); printf("\n");
  printf("d0      = "); printf(NUM_TF, host_a->d0); printf("\n");
  printf("aTau    = "); printf(NUM_TF, host_a->aTau); printf("\n");
  printf("aAlp    = "); printf(NUM_TF, host_a->aAlp); printf("\n");
  printf("aDel    = "); printf(NUM_TF, host_a->aDel); printf("\n");
  printf("bTau    = "); printf(NUM_TF, host_a->bTau); printf("\n");
  printf("bAlp    = "); printf(NUM_TF, host_a->bAlp); printf("\n");
  printf("bDel    = "); printf(NUM_TF, host_a->bDel); printf("\n");
  printf("gamPhi  = "); printf(NUM_TF, host_a->gamPhi); printf("\n");
  printf("gamAlp  = "); printf(NUM_TF, host_a->gamAlp); printf("\n");
  printf("gamDel  = "); printf(NUM_TF, host_a->gamDel); printf("\n");
  printf("sigPhi0 = "); printf(NUM_TF, host_a->sigPhi0); printf("\n");
  printf("sigAlp0 = "); printf(NUM_TF, host_a->sigAlp0); printf("\n");
  printf("sigDel0 = "); printf(NUM_TF, host_a->sigDel0); printf("\n\n");
  
  pf2(host_a->c, host_a->M + 1, host_a->N, "c =\n");
  pf1(host_a->sigC, host_a->M + 1, "sigC =\n");
  pf3(host_a->eps, host_a->M + 1, host_a->N, host_a->G, "eps =\n");
  pf2(host_a->eta, host_a->M + 1, host_a->G, "eta =\n");
  pf1(host_a->d, host_a->M + 1, "d =\n");
  pf1(host_a->tau, host_a->M + 1, "tau =\n");
  pf2(host_a->phi, host_a->M + 1, host_a->G, "phi =\n");
  pf2(host_a->alp, host_a->M + 1, host_a->G, "alp =\n");
  pf2(host_a->del, host_a->M + 1, host_a->G, "del =\n");
  pf1(host_a->thePhi, host_a->M + 1, "thePhi =\n");
  pf1(host_a->theAlp, host_a->M + 1, "theAlp =\n");
  pf1(host_a->theDel, host_a->M + 1, "theDel =\n");
  pf1(host_a->sigPhi, host_a->M + 1, "sigPhi =\n");
  pf1(host_a->sigAlp, host_a->M + 1, "sigAlp =\n");
  pf1(host_a->sigDel, host_a->M + 1, "sigDel =\n");
  pf1(host_a->piAlp, host_a->M + 1, "piAlp =\n");
  pf1(host_a->piDel, host_a->M + 1, "piDel =\n");
  
  printf("mC      = %d\n", host_a->mC);
  printf("mSigC   = %d\n", host_a->mSigC);
  printf("mEps    = %d\n", host_a->mEps);
  printf("mEta    = %d\n", host_a->mEta);
  printf("mD      = %d\n", host_a->mD);
  printf("mTau    = %d\n", host_a->mTau);
  printf("mPhi    = %d\n", host_a->mPhi);
  printf("mAlp    = %d\n", host_a->mAlp);
  printf("mDel    = %d\n", host_a->mDel);
  printf("mThePhi = %d\n", host_a->mThePhi);
  printf("mTheAlp = %d\n", host_a->mTheAlp);
  printf("mTheDel = %d\n", host_a->mTheDel);
  printf("mSigPhi = %d\n", host_a->mSigPhi);
  printf("mSigAlp = %d\n", host_a->mSigAlp);
  printf("mSigDel = %d\n", host_a->mSigDel);
  printf("mPiAlp  = %d\n", host_a->mPiAlp);
  printf("mPiDel  = %d\n\n", host_a->mPiDel);
  
  printf("tuneD = "); printf(NUM_TF, host_a->tuneD); printf("\n");
  pf1(host_a->tuneC, host_a->N, "tuneC =\n");
  pf1(host_a->tunePhi, host_a->G, "tunePhi =\n");
  pf2(host_a->tuneEps, host_a->N, host_a->G, "tuneEps =\n");
  
  printf("accD = %d\n", host_a->accD);
  pi1(host_a->accC, host_a->N, "accC =\n");
  pi1(host_a->accPhi, host_a->G, "accPhi =\n");
  pi1(host_a->accAlp, host_a->G, "accAlp =\n");
  pi1(host_a->accDel, host_a->G, "accDel =\n");
  pi2(host_a->accEps, host_a->N, host_a->G, "accEps =\n");
    
  printf("s1 = "); printf(NUM_TF, host_a->s1); printf("\n");
  printf("s2 = "); printf(NUM_TF, host_a->s2); printf("\n\n");
  
  pf1(host_a->tmp1, host_a->G, "tmp1 =\n");
  pf1(host_a->tmp1, host_a->G, "tmp2 =\n");
  
  pf1(host_a->Old, host_a->N, "Old =\n");
  pf1(host_a->New, host_a->N, "New =\n");
  pf1(host_a->lOld, host_a->N, "lOld =\n");
  pf1(host_a->lNew, host_a->N, "lNew =\n");
  
  printf("constSigC = %d\n", host_a->constSigC);
  printf("constD = %d\n", host_a->constD);
  printf("constTau = %d\n", host_a->constTau);
  printf("constThePhi = %d\n", host_a->constThePhi);
  printf("constTheAlp = %d\n", host_a->constTheAlp);
  printf("constTheDel = %d\n", host_a->constTheDel);
  printf("constSigPhi = %d\n", host_a->constSigPhi);
  printf("constSigAlp = %d\n", host_a->constSigAlp);
  printf("constSigDel = %d\n", host_a->constSigDel);
  printf("constPiAlp = %d\n", host_a->constPiAlp);
  printf("constPiDel = %d\n", host_a->constPiDel);
  
  if(!onHost)
    freeChain(host_a, cfg, 1);
}