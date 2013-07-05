#include <Config.h>
#include <stdio.h>
#include <stdlib.h>

void printConfig(Config *cfg){
  
  printf("dataFile = %s\n", cfg->dataFile);
  printf("groupFile = %s\n", cfg->dataGroup);
  printf("probsFile = %s\n", cfg->probsFile);
  printf("hyperFile = %s\n", cfg->hyperFile);
  printf("ratesFile = %s\n", cfg->ratesFile);
  printf("parmsFile = %s\n\n", cfg->parmsFile);

  printf("probsFlag = %d\n", cfg->probsFlag);
  printf("hyperFlag = %d\n", cfg->hyperFlag);
  printf("ratesFlag = %d\n", cfg->ratesFlag);
  printf("parmsFlag = %d\n\n", cfg->parmsFlag);
  
  printf("burnin = %d\n", cfg->burnin);
  printf("joint = %d\n", cfg->joint);
  printf("heterosis = %d\n\n", cfg->heterosis);
  
  printf("M = %d\n", cfg->M);
  printf("N = %d\n", cfg->N);
  printf("G = %d\n\n", cfg->G);
  
  printf("sigC0 = %0.3f\n", cfg->sigC0);
  printf("d0 = %0.3f\n", cfg->d0);
  printf("aTau = %0.3f\n", cfg->aTau);
  printf("aAlp = %0.3f\n", cfg->aAlp);
  printf("aDel = %0.3f\n", cfg->aDel);
  printf("bTau = %0.3f\n", cfg->bTau);
  printf("bAlp = %0.3f\n", cfg->bAlp);
  printf("bDel = %0.3f\n", cfg->bDel);
  printf("gamPhi = %0.3f\n", cfg->gamPhi);
  printf("gamAlp = %0.3f\n", cfg->gamAlp);
  printf("gamDel = %0.3f\n", cfg->gamDel);
  printf("sigPhi0 = %0.3f\n", cfg->sigPhi0);
  printf("sigAlp0 = %0.3f\n", cfg->sigAlp0);
  printf("sigDel0 = %0.3f\n\n", cfg->sigDel0);
  
  printf("sigC = %0.3f\n", cfg->sigC);
  printf("d = %0.3f\n", cfg->d);
  printf("tau = %0.3f\n", cfg->tau);
  printf("thePhi = %0.3f\n", cfg->thePhi);
  printf("theAlp = %0.3f\n", cfg->theAlp);
  printf("theDel = %0.3f\n", cfg->theDel);
  printf("sigPhi = %0.3f\n", cfg->sigPhi);
  printf("sigAlp = %0.3f\n", cfg->sigAlp);
  printf("sigDel = %0.3f\n", cfg->sigDel);
  printf("piAlp = %0.3f\n", cfg->piAlp);
  printf("piDel = %0.3f\n\n", cfg->piDel);
  
  printf("constD = %0.3f\n", cfg->constD);
  printf("constTau = %0.3f\n", cfg->constTau);
  printf("constThePhi = %0.3f\n", cfg->constThePhi);
  printf("constTheAlp = %0.3f\n", cfg->constTheAlp);
  printf("constTheDel = %0.3f\n", cfg->constTheDel);
  printf("constSigPhi = %0.3f\n", cfg->constSigPhi);
  printf("constSigAlp = %0.3f\n", cfg->constSigAlp);
  printf("constSigDel = %0.3f\n", cfg->constSigDel);
  printf("constPiAlp = %0.3f\n", cfg->constPiAlp);
  printf("constPiDel = %0.3f\n", cfg->constPiDel);
}