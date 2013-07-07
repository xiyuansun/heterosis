#include <Chain.h>
#include <Config.h>
#include <functions.h>
#include <numericTypes.h>
#include <stdio.h>
#include <stdlib.h>

void printChain(Chain *a){

  printf("M = %d\n", a->M);
  printf("N = %d\n", a->N);
  printf("G = %d\n\n", a->G);  
  
  printf("burnin = %d\n", a->burnin);
  printf("heterosis = %d\n", a->heterosis);
  printf("someParmsFlag = %d\n", a->someParmsFlag);
  printf("allParmsFlag = %d\n\n", a->allParmsFlag); 

  pi2(a->y, a->N, a->G, "y = \n");
  pf1(a->yMeanG, a->N, "yMeanG =\n");
  pi1(a->grp, a->N, "grp =\n");
  
  printf("sigC0   = %0.3f\n", a->sigC0);
  printf("d0      = %0.3f\n", a->d0);
  printf("aTau    = %0.3f\n", a->aTau);
  printf("aAlp    = %0.3f\n", a->aAlp);
  printf("aDel    = %0.3f\n", a->aDel);
  printf("bTau    = %0.3f\n", a->bTau);
  printf("bAlp    = %0.3f\n", a->bAlp);
  printf("bDel    = %0.3f\n", a->bDel);  
  printf("gamPhi  = %0.3f\n", a->gamPhi);
  printf("gamAlp  = %0.3f\n", a->gamAlp);
  printf("gamDel  = %0.3f\n", a->gamDel);
  printf("sigPhi0 = %0.3f\n", a->sigPhi0);
  printf("sigAlp0 = %0.3f\n", a->sigAlp0);
  printf("sigDel0 = %0.3f\n\n", a->sigDel0);
  
  pf2(a->c, a->M, a->N, "c =\n");
  pf1(a->sigC, a->M, "sigC =\n");
  pf3(a->eps, a->M, a->N, a->G, "eps =\n");
  pf2(a->eta, a->M, a->N, "eta =\n");
  pf1(a->d, a->M, "d =\n");
  pf1(a->tau, a->M, "tau =\n");
  pf2(a->phi, a->M, a->G, "phi =\n");
  pf2(a->alp, a->M, a->G, "alp =\n");
  pf2(a->del, a->M, a->G, "del =\n");
  pf1(a->thePhi, a->M, "thePhi =\n");
  pf1(a->theAlp, a->M, "theAlp =\n");
  pf1(a->theDel, a->M, "theDel =\n");
  pf1(a->sigPhi, a->M, "sigPhi =\n");
  pf1(a->sigAlp, a->M, "sigAlp =\n");
  pf1(a->sigDel, a->M, "sigDel =\n");
  pf1(a->piAlp, a->M, "piAlp =\n");
  pf1(a->piDel, a->M, "piDel =\n");
  
  printf("mC      = %d\n", a->mC);
  printf("mSigC   = %d\n", a->mSigC);
  printf("mEta    = %d\n", a->mEta);
  printf("mD      = %d\n", a->mD);
  printf("mTau    = %d\n", a->mTau);
  printf("mPhi    = %d\n", a->mPhi);
  printf("mThePhi = %d\n", a->mThePhi);
  printf("mSigPhi = %d\n", a->mSigPhi);
  printf("mAlp    = %d\n", a->mAlp);
  printf("mTheAlp = %d\n", a->mTheAlp);
  printf("mSigAlp = %d\n", a->mSigAlp);
  printf("mPiAlp  = %d\n", a->mPiAlp);
  printf("mDel    = %d\n", a->mDel);
  printf("mTheDel = %d\n", a->mTheDel);
  printf("mSigDel = %d\n", a->mSigDel);
  printf("mPiDel  = %d\n\n", a->mPiDel);
  
  printf("tuneD = %0.3f\n", a->tuneD);
  pf1(a->tuneC, a->N, "tuneC =\n");
  pf1(a->tunePhi, a->G, "tunePhi =\n");
  pf2(a->tuneEps, a->N, a->G, "tuneEps =\n");
  
  printf("accD = %d\n", a->accD);
  pi1(a->accC, a->N, "accC =\n");
  pi1(a->accPhi, a->G, "accPhi =\n");
  pi1(a->accAlp, a->G, "accAlp =\n");
  pi1(a->accDel, a->G, "accDel =\n");
  pi2(a->accEps, a->N, a->G, "accEps =\n");
    
  printf("shape = %0.3f\n", a->shape);
  printf("rate = %0.3f\n\n", a->rate);
  
  pf1(a->tmp1, a->G, "tmp1 =\n");
  pf1(a->tmp1, a->G, "tmp2 =\n");
  
  pf2(a->Old, a->N, a->G, "Old =\n");
  pf2(a->New, a->N, a->G, "New =\n");
  pf2(a->lOld, a->N, a->G, "lOld =\n");
  pf2(a->lNew, a->N, a->G, "lNew =\n");
}