#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

__host__ void printChain(Chain *host_a, Chain *dev_a, Config *cfg){
  int N = cfg->N, G = cfg->G;
  Chain *allHost_a = chainDeviceToHost(host_a, dev_a, cfg);

  printf("M = %d\n", allHost_a->M); ;
  printf("N = %d\n", allHost_a->N);
  printf("G = %d\n\n", allHost_a->G);  
  
  printf("burnin = %d\n", allHost_a->burnin);
  printf("heterosis = %d\n", allHost_a->heterosis);
  printf("parmsFlag = %d\n", allHost_a->parmsFlag);

  pi2(allHost_a->y, allHost_a->N, allHost_a->G, "y = \n");
  pf1(allHost_a->yMeanG, allHost_a->N, "yMeanG =\n");
  pi1(allHost_a->grp, allHost_a->N, "grp =\n");  
  pstate(allHost_a->states, MAX_NG, "curand states =\n");
  
  printf("sigC0   = "); printf(NUM_TF, allHost_a->sigC0); printf("\n");
  printf("d0      = "); printf(NUM_TF, allHost_a->d0); printf("\n");
  printf("aTau    = "); printf(NUM_TF, allHost_a->aTau); printf("\n");
  printf("aAlp    = "); printf(NUM_TF, allHost_a->aAlp); printf("\n");
  printf("aDel    = "); printf(NUM_TF, allHost_a->aDel); printf("\n");
  printf("bTau    = "); printf(NUM_TF, allHost_a->bTau); printf("\n");
  printf("bAlp    = "); printf(NUM_TF, allHost_a->bAlp); printf("\n");
  printf("bDel    = "); printf(NUM_TF, allHost_a->bDel); printf("\n");
  printf("gamPhi  = "); printf(NUM_TF, allHost_a->gamPhi); printf("\n");
  printf("gamAlp  = "); printf(NUM_TF, allHost_a->gamAlp); printf("\n");
  printf("gamDel  = "); printf(NUM_TF, allHost_a->gamDel); printf("\n");
  printf("sigPhi0 = "); printf(NUM_TF, allHost_a->sigPhi0); printf("\n");
  printf("sigAlp0 = "); printf(NUM_TF, allHost_a->sigAlp0); printf("\n");
  printf("sigDel0 = "); printf(NUM_TF, allHost_a->sigDel0); printf("\n\n");
  
  pf2(allHost_a->c, allHost_a->M + 1, allHost_a->N, "c =\n");
  pf1(allHost_a->sigC, allHost_a->M + 1, "sigC =\n");
  pf3(allHost_a->eps, allHost_a->M + 1, allHost_a->N, allHost_a->G, "eps =\n");
  pf2(allHost_a->eta, allHost_a->M + 1, allHost_a->G, "eta =\n");
  pf1(allHost_a->d, allHost_a->M + 1, "d =\n");
  pf1(allHost_a->tau, allHost_a->M + 1, "tau =\n");
  pf2(allHost_a->phi, allHost_a->M + 1, allHost_a->G, "phi =\n");
  pf2(allHost_a->alp, allHost_a->M + 1, allHost_a->G, "alp =\n");
  pf2(allHost_a->del, allHost_a->M + 1, allHost_a->G, "del =\n");
  pf1(allHost_a->thePhi, allHost_a->M + 1, "thePhi =\n");
  pf1(allHost_a->theAlp, allHost_a->M + 1, "theAlp =\n");
  pf1(allHost_a->theDel, allHost_a->M + 1, "theDel =\n");
  pf1(allHost_a->sigPhi, allHost_a->M + 1, "sigPhi =\n");
  pf1(allHost_a->sigAlp, allHost_a->M + 1, "sigAlp =\n");
  pf1(allHost_a->sigDel, allHost_a->M + 1, "sigDel =\n");
  pf1(allHost_a->piAlp, allHost_a->M + 1, "piAlp =\n");
  pf1(allHost_a->piDel, allHost_a->M + 1, "piDel =\n");
  
  printf("mC      = %d\n", allHost_a->mC);
  printf("mSigC   = %d\n", allHost_a->mSigC);
  printf("mEps    = %d\n", allHost_a->mEps);
  printf("mEta    = %d\n", allHost_a->mEta);
  printf("mD      = %d\n", allHost_a->mD);
  printf("mTau    = %d\n", allHost_a->mTau);
  printf("mPhi    = %d\n", allHost_a->mPhi);
  printf("mAlp    = %d\n", allHost_a->mAlp);
  printf("mDel    = %d\n", allHost_a->mDel);
  printf("mThePhi = %d\n", allHost_a->mThePhi);
  printf("mTheAlp = %d\n", allHost_a->mTheAlp);
  printf("mTheDel = %d\n", allHost_a->mTheDel);
  printf("mSigPhi = %d\n", allHost_a->mSigPhi);
  printf("mSigAlp = %d\n", allHost_a->mSigAlp);
  printf("mSigDel = %d\n", allHost_a->mSigDel);
  printf("mPiAlp  = %d\n", allHost_a->mPiAlp);
  printf("mPiDel  = %d\n\n", allHost_a->mPiDel);
  
  printf("tuneD = "); printf(NUM_TF, allHost_a->tuneD); printf("\n");
  pf1(allHost_a->tuneC, allHost_a->N, "tuneC =\n");
  pf1(allHost_a->tunePhi, allHost_a->G, "tunePhi =\n");
  pf2(allHost_a->tuneEps, allHost_a->N, allHost_a->G, "tuneEps =\n");
  
  printf("accD = %d\n", allHost_a->accD);
  pi1(allHost_a->accC, allHost_a->N, "accC =\n");
  pi1(allHost_a->accPhi, allHost_a->G, "accPhi =\n");
  pi1(allHost_a->accAlp, allHost_a->G, "accAlp =\n");
  pi1(allHost_a->accDel, allHost_a->G, "accDel =\n");
  pi2(allHost_a->accEps, allHost_a->N, allHost_a->G, "accEps =\n");
    
  printf("s1 = "); printf(NUM_TF, allHost_a->s1); printf("\n");
  printf("s2 = "); printf(NUM_TF, allHost_a->s2); printf("\n\n");
  
  pf1(allHost_a->tmp1, allHost_a->G, "tmp1 =\n");
  pf1(allHost_a->tmp1, allHost_a->G, "tmp2 =\n");
  
  pf1(allHost_a->Old, allHost_a->N, "Old =\n");
  pf1(allHost_a->New, allHost_a->N, "New =\n");
  pf1(allHost_a->lOld, allHost_a->N, "lOld =\n");
  pf1(allHost_a->lNew, allHost_a->N, "lNew =\n");
  
  printf("constSigC = %d\n", allHost_a->constSigC);
  printf("constD = %d\n", allHost_a->constD);
  printf("constTau = %d\n", allHost_a->constTau);
  printf("constThePhi = %d\n", allHost_a->constThePhi);
  printf("constTheAlp = %d\n", allHost_a->constTheAlp);
  printf("constTheDel = %d\n", allHost_a->constTheDel);
  printf("constSigPhi = %d\n", allHost_a->constSigPhi);
  printf("constSigAlp = %d\n", allHost_a->constSigAlp);
  printf("constSigDel = %d\n", allHost_a->constSigDel);
  printf("constPiAlp = %d\n", allHost_a->constPiAlp);
  printf("constPiDel = %d\n", allHost_a->constPiDel);
  
  freeChain(allHost_a, cfg, 1);
}