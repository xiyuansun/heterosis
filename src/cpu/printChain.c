#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

void printChain(Chain *a){
  printf("Chain object:\n\n");

  printf("m = %d\n", a->m);
  printf("M = %d\n", a->M);
  printf("N = %d\n", a->N);
  printf("G = %d\n\n", a->G);  
  printf("burnin = %d\n", a->burnin);

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
  
  pf1(a->c, a->N, "c =\n");
  printf("sigC = "); printf(NUM_TF, a->sigC); printf("\n\n"); 
  pf2(a->eps, a->N, a->G, "eps =\n");
  pf1(a->eta, a->G, "eta =\n");
  printf("d = "); printf(NUM_TF, a->d); printf("\n");
  printf("tau = "); printf(NUM_TF, a->tau); printf("\n");
  pf1(a->phi, a->G, "phi =\n");
  pf1(a->alp, a->G, "alp =\n");
  pf1(a->del, a->G, "del =\n");
  printf("thePhi = "); printf(NUM_TF, a->thePhi); printf("\n");
  printf("theAlp = "); printf(NUM_TF, a->theAlp); printf("\n");
  printf("theDel = "); printf(NUM_TF, a->theDel); printf("\n");
  printf("sigPhi = "); printf(NUM_TF, a->sigPhi); printf("\n");
  printf("sigAlp = "); printf(NUM_TF, a->sigAlp); printf("\n");
  printf("sigDel = "); printf(NUM_TF, a->sigDel); printf("\n");
  printf("piAlp = "); printf(NUM_TF, a->piAlp); printf("\n"); 
  printf("piDel = "); printf(NUM_TF, a->piDel); printf("\n\n");
  
  pi1(a->dex, a->G, "dex =\n");
  pi1(a->hph, a->G, "hph =\n");
  pi1(a->lph, a->G, "lph =\n"); 
  pi1(a->mph, a->G, "mph =\n");
  
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
    
  printf("s1 = "); printf(NUM_TF, a->s1); printf("\n");
  printf("s2 = "); printf(NUM_TF, a->s2); printf("\n\n");
  
  pf1(a->tmp1, a->G, "tmp1 =\n");
  pf1(a->tmp1, a->G, "tmp2 =\n");
  
  pf1(a->Old, a->N, "Old =\n");
  pf1(a->New, a->N, "New =\n");
  pf1(a->lOld, a->N, "lOld =\n");
  pf1(a->lNew, a->N, "lNew =\n");
  
  printf("constSigC = %d\n", a->constSigC);
  printf("constD = %d\n", a->constD);
  printf("constTau = %d\n", a->constTau);
  printf("constThePhi = %d\n", a->constThePhi);
  printf("constTheAlp = %d\n", a->constTheAlp);
  printf("constTheDel = %d\n", a->constTheDel);
  printf("constSigPhi = %d\n", a->constSigPhi);
  printf("constSigAlp = %d\n", a->constSigAlp);
  printf("constSigDel = %d\n", a->constSigDel);
  printf("constPiAlp = %d\n", a->constPiAlp);
  printf("constPiDel = %d\n\n", a->constPiDel);

  /* for computing DIC */
  
  printf("meanLogLik = "); printf(NUM_TF, a->meanLogLik); printf("\n"); 
  printf("logLikMeans = "); printf(NUM_TF, a->logLikMean); printf("\n"); 
  printf("dic = "); printf(NUM_TF, a->dic); printf("\n"); 
  
  pf1(a->meanC, a->N, "meanC =\n");
  pf1(a->meanPhi, a->G, "meanPhi =\n");
  pf1(a->meanAlp, a->G, "meanAlp =\n");
  pf1(a->meanDel, a->G, "meanDel =\n"); 
  pf2(a->meanEps, a->N, a->G, "meanEps =\n");

}