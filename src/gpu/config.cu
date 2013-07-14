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
  cfg->timeFlag = 0;
  cfg->verbose = 0;
  cfg->diagnostics = 0;

  cfg->chains = 2;
  cfg->m = 1;
  cfg->M = 10;
  cfg->N = 0;
  cfg->G = 0;
  cfg->burnin = cfg->M / 2;
  cfg->joint = 0;
  cfg->seed = 22;

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
  
  cfg->timeC = 0;
  cfg->timeTau = 0;
  cfg->timePiAlp = 0;
  cfg->timePiDel = 0;
  cfg->timeD = 0;
  cfg->timeThePhi = 0;
  cfg->timeTheAlp = 0;
  cfg->timeTheDel = 0;
  cfg->timeSigC = 0;
  cfg->timeSigPhi = 0;
  cfg->timeSigAlp = 0;
  cfg->timeSigDel = 0;
  cfg->timeEta = 0;
  cfg->timeEps = 0;
  cfg->timePhi = 0;
  cfg->timeAlp = 0;
  cfg->timeDel = 0;

  getopts(cfg, argc, argv);
  srand(cfg->seed);
   
  system("mkdir -p ../out/");
  system("mkdir -p ../out/probs/");
  
  if(cfg->ratesFlag)
    system("mkdir -p ../out/rates/");
  
  if(cfg->hyperFlag)
    system("mkdir -p ../out/hyper/");
  
  if(cfg->parmsFlag)
    system("mkdir -p ../out/parms/"); 
  
  if(cfg->timeFlag)
    system("mkdir -p ../out/time/");
  
  if(cfg->diagnostics){
    system("mkdir -p ../out/diagnostics/");
    system("rm -f ../out/diagnostics/dic.txt");
  }
  
  return cfg;
}