#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
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
  
  printf("sigC0   = "); printf(NUM_TF, a->sigC0); printf("\n");
  printf("d0      = "); printf(NUM_TF, a->d0); printf("\n");
  printf("aTau    = "); printf(NUM_TF, a->aTau); printf("\n");
  printf("aAlp    = "); printf(NUM_TF, a->aAlp); printf("\n");
  printf("aDel    = "); printf(NUM_TF, a->aDel); printf("\n");
  printf("bTau    = "); printf(NUM_TF, a->bTau); printf("\n");
  printf("bAlp    = "); printf(NUM_TF, a->bAlp); printf("\n");
  printf("bDel    = "); printf(NUM_TF, a->bDel); printf("\n");
  printf("gamPhi  = "); printf(NUM_TF, a->gamPhi); printf("\n");
  printf("gamAlp  = "); printf(NUM_TF, a->gamAlp); printf("\n");
  printf("gamDel  = "); printf(NUM_TF, a->gamDel); printf("\n");
  printf("sigPhi0 = "); printf(NUM_TF, a->sigPhi0); printf("\n");
  printf("sigAlp0 = "); printf(NUM_TF, a->sigAlp0); printf("\n");
  printf("sigDel0 = "); printf(NUM_TF, a->sigDel0); printf("\n\n");
  
  pf2(a->c, a->M + 1, a->N, "c =\n");
  pf1(a->sigC, a->M + 1, "sigC =\n");
  pf3(a->eps, a->M + 1, a->N, a->G, "eps =\n");
  pf2(a->eta, a->M + 1, a->N, "eta =\n");
  pf1(a->d, a->M + 1, "d =\n");
  pf1(a->tau, a->M + 1, "tau =\n");
  pf2(a->phi, a->M + 1, a->G, "phi =\n");
  pf2(a->alp, a->M + 1, a->G, "alp =\n");
  pf2(a->del, a->M + 1, a->G, "del =\n");
  pf1(a->thePhi, a->M + 1, "thePhi =\n");
  pf1(a->theAlp, a->M + 1, "theAlp =\n");
  pf1(a->theDel, a->M + 1, "theDel =\n");
  pf1(a->sigPhi, a->M + 1, "sigPhi =\n");
  pf1(a->sigAlp, a->M + 1, "sigAlp =\n");
  pf1(a->sigDel, a->M + 1, "sigDel =\n");
  pf1(a->piAlp, a->M + 1, "piAlp =\n");
  pf1(a->piDel, a->M + 1, "piDel =\n");
  
  printf("mC      = %d\n", a->mC);
  printf("mSigC   = %d\n", a->mSigC);
  printf("mEps    = %d\n", a->mEps);
  printf("mEta    = %d\n", a->mEta);
  printf("mD      = %d\n", a->mD);
  printf("mTau    = %d\n", a->mTau);
  printf("mPhi    = %d\n", a->mPhi);
  printf("mAlp    = %d\n", a->mAlp);
  printf("mDel    = %d\n", a->mDel);
  printf("mThePhi = %d\n", a->mThePhi);
  printf("mTheAlp = %d\n", a->mTheAlp);
  printf("mTheDel = %d\n", a->mTheDel);
  printf("mSigPhi = %d\n", a->mSigPhi);
  printf("mSigAlp = %d\n", a->mSigAlp);
  printf("mSigDel = %d\n", a->mSigDel);
  printf("mPiAlp  = %d\n", a->mPiAlp);
  printf("mPiDel  = %d\n\n", a->mPiDel);
  
  printf("tuneD = "); printf(NUM_TF, a->tuneD); printf("\n");
  pf1(a->tuneC, a->N, "tuneC =\n");
  pf1(a->tunePhi, a->G, "tunePhi =\n");
  pf2(a->tuneEps, a->N, a->G, "tuneEps =\n");
  
  printf("accD = %d\n", a->accD);
  pi1(a->accC, a->N, "accC =\n");
  pi1(a->accPhi, a->G, "accPhi =\n");
  pi1(a->accAlp, a->G, "accAlp =\n");
  pi1(a->accDel, a->G, "accDel =\n");
  pi2(a->accEps, a->N, a->G, "accEps =\n");
    
  printf("shape = "); printf(NUM_TF, a->shape); printf("\n");
  printf("rate = "); printf(NUM_TF, a->rate); printf("\n\n");
    
  printf("s1 = "); printf(NUM_TF, a->s1); printf("\n");
  printf("s2 = "); printf(NUM_TF, a->s2); printf("\n\n");
  
  pf1(a->tmp1, a->G, "tmp1 =\n");
  pf1(a->tmp1, a->G, "tmp2 =\n");
  
  pf2(a->Old, a->N, a->G, "Old =\n");
  pf2(a->New, a->N, a->G, "New =\n");
  pf2(a->lOld, a->N, a->G, "lOld =\n");
  pf2(a->lNew, a->N, a->G, "lNew =\n");
}