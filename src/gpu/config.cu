#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Config *config(int argc, char **argv){

  Config *cfg = (Config*) malloc(sizeof(Config));
  cfg->chainNum = 1;
  
  /* default filenames */        

  strcpy(cfg->dataFile, "../data/data.txt"); 
  strcpy(cfg->groupFile, "../data/group.txt");
   
  cfg->ratesFlag = 0;
  cfg->hyperFlag = 0;
  cfg->parmsFlag = 0;

  cfg->M = 10;  
  cfg->burnin = cfg->M / 2;
  cfg->joint = 0;

  /* default initialization constants */

  cfg->sigC0 = 10;
  cfg->d0 = 1000;
  cfg->aTau = 100;
  cfg->aAlp = 1;
  cfg->aDel = 1;
  cfg->bTau = 100;
  cfg->bAlp = 1;
  cfg->bDel = 1;
  cfg->gamPhi = 2;
  cfg->gamAlp = 2;
  cfg->gamDel = 2;
  cfg->sigPhi0 = 2;
  cfg->sigAlp0 = 2;
  cfg->sigDel0 = 2;
  
  /* default: don't hold hyperparameters constant */
  
  cfg->constSigC = 0;
  cfg->constD = 0;
  cfg->constTau = 0;
  cfg->constThePhi = 0;
  cfg->constTheAlp = 0;
  cfg->constTheDel = 0;
  cfg->constSigPhi = 0;
  cfg->constSigAlp = 0;
  cfg->constSigDel = 0;
  cfg->constPiAlp = 0;
  cfg->constPiDel = 0;

  getopts(cfg, argc, argv);
  srand(cfg->seed);

  /* initialize output and chain */

  system("mkdir -p ../out/");
  system("mkdir -p ../out/probs/");
  
  cfg->log = fopen("../out/log.txt", "w+");
  
  cfg->time = fopen("../out/time.txt", "w+");
  cfg->timeConfig = fopen("../out/timeConfig.txt", "w+");
  fprintf(cfg->time, "c ");

  if(!cfg->constTau){
    cfg->tau = sqrt(rgamma(cfg->aTau, cfg->bTau, 0));
    fprintf(cfg->time, "tau ");
  }
 
  if(!cfg->constPiAlp){
    cfg->piAlp = rbeta(cfg->aAlp, cfg->bAlp);
    fprintf(cfg->time, "piDel ");
  }
  
  if(!cfg->constPiDel){
    cfg->piDel = rbeta(cfg->aDel, cfg->bDel);
    fprintf(cfg->time, "piDel ");
  }
 
  if(!cfg->constD){
    cfg->d = runiform(0, cfg->d0);
    fprintf(cfg->time, "d ");
  }
 
  if(!cfg->constThePhi){
    cfg->thePhi = rnormal(0, cfg->gamPhi);
    fprintf(cfg->time, "thePhi ");
  }

  if(!cfg->constTheAlp){
    cfg->theAlp = rnormal(0, cfg->gamAlp);
    fprintf(cfg->time, "theAlp ");
  }

  if(!cfg->constTheDel){
    cfg->theDel = rnormal(0, cfg->gamDel);
    fprintf(cfg->time, "theDel ");
  }
 
  if(!cfg->constSigC){
    cfg->sigC = runiform(0, cfg->sigC0);
    fprintf(cfg->time, "sigC ");
  }
  
  if(!cfg->constSigPhi){
    cfg->sigPhi = runiform(0, cfg->sigPhi0);
    fprintf(cfg->time, "sigPhi ");
  }

  if(!cfg->constSigAlp){
    cfg->sigAlp = runiform(0, cfg->sigAlp0);
    fprintf(cfg->time, "sigAlp ");
  }

  if(!cfg->constSigDel){
    cfg->sigDel = runiform(0, cfg->sigDel0);
    fprintf(cfg->time, "sigDel ");
  }
  
  fprintf(cfg->time, "eta eps ");
  
  if(cfg->joint){
    fprintf(cfg->time, "phiAlpDel ");
  } else {
    fprintf(cfg->time, "phi alp del ");
  }
  
  fprintf(cfg->timeConfig, "alloc probs ");
  
  if(cfg->ratesFlag){
    system("mkdir -p ../out/rates/");
    fprintf(cfg->timeConfig, "rates ");
  }
  
  if(cfg->hyperFlag){
    system("mkdir -p ../out/hyper/");
    fprintf(cfg->timeConfig, "hyper ");
  }
  
  if(cfg->parmsFlag){
    system("mkdir -p ../out/parms/"); 
    fprintf(cfg->timeConfig, "parms ");
  }
  
  fprintf(cfg->timeConfig, "\n");
  fprintf(cfg->time, "\n"); 
  return cfg;
}