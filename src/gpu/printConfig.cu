#include <Config.h>
#include <constants.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void printConfig(Config *cfg){

  printf("chainNum = %d\n", cfg->chainNum);
  printf("m = %d\n", cfg->m);
  printf("chains = %d\n\n", cfg->chains);
  
  printf("dataFile = %s\n", cfg->dataFile);
  printf("groupFile = %s\n", cfg->groupFile);

  printf("probs = %d\n", cfg->probs);
  printf("rates = %d\n", cfg->rates);
  printf("hyper = %d\n", cfg->hyper);
  printf("parms = %d\n", cfg->parms);
  printf("time = %d\n", cfg->time);
  printf("verbose = %d\n", cfg->verbose);
  printf("dic = %d\n\n", cfg->dic);
 
  printf("burnin = %d\n", cfg->burnin);
  printf("joint = %d\n", cfg->joint);
  printf("seed = %d\n\n", cfg->seed);
  printf("heterosis = %d\n\n", cfg->heterosis);
  
  printf("M = %d\n", cfg->M);
  printf("N = %d\n", cfg->N);
  printf("G = %d\n\n", cfg->G);
  
  printf("sigC0 = "); printf(NUM_TF, cfg->sigC0); printf("\n");
  printf("d0 = "); printf(NUM_TF, cfg->d0); printf("\n");
  printf("aTau = "); printf(NUM_TF, cfg->aTau); printf("\n");
  printf("aAlp = "); printf(NUM_TF, cfg->aAlp); printf("\n"); 
  printf("aDel = "); printf(NUM_TF, cfg->aDel); printf("\n");
  printf("bTau = "); printf(NUM_TF, cfg->bTau); printf("\n");
  printf("bAlp = "); printf(NUM_TF, cfg->bAlp); printf("\n");
  printf("bDel = "); printf(NUM_TF, cfg->bDel); printf("\n");
  printf("gamPhi = "); printf(NUM_TF, cfg->gamPhi); printf("\n");
  printf("gamAlp = "); printf(NUM_TF, cfg->gamAlp); printf("\n");
  printf("gamDel = "); printf(NUM_TF, cfg->gamDel); printf("\n");
  printf("sigPhi0 = "); printf(NUM_TF, cfg->sigPhi0); printf("\n");
  printf("sigAlp0 = "); printf(NUM_TF, cfg->sigAlp0); printf("\n");
  printf("sigDel0 = "); printf(NUM_TF, cfg->sigDel0); printf("\n\n");
  
  printf("sigC = "); printf(NUM_TF, cfg->sigC); printf("\n");
  printf("d = "); printf(NUM_TF, cfg->d); printf("\n");
  printf("tau = "); printf(NUM_TF, cfg->tau); printf("\n");
  printf("thePhi = "); printf(NUM_TF, cfg->thePhi); printf("\n");
  printf("theAlp = "); printf(NUM_TF, cfg->theAlp); printf("\n");
  printf("theDel = "); printf(NUM_TF, cfg->theDel); printf("\n");
  printf("sigPhi = "); printf(NUM_TF, cfg->sigPhi); printf("\n");
  printf("sigAlp = "); printf(NUM_TF, cfg->sigAlp); printf("\n");
  printf("sigDel = "); printf(NUM_TF, cfg->sigDel); printf("\n");
  printf("piAlp = "); printf(NUM_TF, cfg->piAlp); printf("\n");
  printf("piDel = "); printf(NUM_TF, cfg->piDel); printf("\n\n");
  
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