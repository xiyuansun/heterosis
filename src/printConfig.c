#include <Config.h>
#include <constants.h>
#include <stdio.h>
#include <stdlib.h>

void printConfig(Config *cfg){
  
  printf("dataFile = %s\n", cfg->dataFile);
  printf("groupFile = %s\n", cfg->groupFile);
  printf("probsFile = %s\n", cfg->probsFile);
  printf("ratesFile = %s\n", cfg->ratesFile);
  printf("hyperFile = %s\n", cfg->hyperFile);
  printf("someParmsFile = %s\n", cfg->someParmsFile);
  printf("allParmsFile = %s\n\n", cfg->allParmsFile);

  printf("probsFlag = %d\n", cfg->probsFlag);
  printf("ratesFlag = %d\n", cfg->ratesFlag);
  printf("hyperFlag = %d\n", cfg->hyperFlag);
  printf("someParmsFlag = %d\n", cfg->someParmsFlag);
  printf("allParmsFlag = %d\n\n", cfg->allParmsFlag);  
  
  printf("burnin = %d\n", cfg->burnin);
  printf("joint = %d\n", cfg->joint);
  printf("seed = %d\n\n", cfg->seed);
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
  
  printf("constSigC = %d\n", cfg->constSigC);
  printf("constD = %d\n", cfg->constD);
  printf("constTau = %d\n", cfg->constTau);
  printf("constThePhi = %d\n", cfg->constThePhi);
  printf("constTheAlp = %d\n", cfg->constTheAlp);
  printf("constTheDel = %d\n", cfg->constTheDel);
  printf("constSigPhi = %d\n", cfg->constSigPhi);
  printf("constSigAlp = %d\n", cfg->constSigAlp);
  printf("constSigDel = %d\n", cfg->constSigDel);
  printf("constPiAlp = %d\n", cfg->constPiAlp);
  printf("constPiDel = %d\n", cfg->constPiDel);
}