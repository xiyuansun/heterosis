#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

Config *config(int argc, char **argv){

  Config *cfg = (Config*) malloc(sizeof(Config));
  mode_t process_mask;
  
  /* default filenames */        

  strcpy(cfg->dataFile, "../data/data.txt"); 
  strcpy(cfg->groupFile, "../data/group.txt");
  strcpy(cfg->outDir, "../out");
   
  cfg->probs = 0; 
  cfg->rates = 0;
  cfg->hyper = 0;
  cfg->parms = 0;
  cfg->time = 0;
  cfg->verbose = 0;
  cfg->dic = 0;

  cfg->chains = 2;
  cfg->chainNum = 0;
  cfg->m = 1;
  cfg->M = 10;
  cfg->N = 0;
  cfg->G = 0;
  cfg->burnin = cfg->M / 2;
  cfg->joint = 0;
  cfg->seed = 22;
  cfg->heterosis = 1;

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
  
  /* unset hyperparameters get value -1 */
  
  cfg->sigC = -1;
  cfg->d = -1;
  cfg->tau = -1;
  cfg->thePhi = -1;
  cfg->theAlp = -1;
  cfg->theDel = -1;
  cfg->sigPhi = -1;
  cfg->sigAlp = -1;
  cfg->sigDel = -1;
  cfg->piAlp = -1;
  cfg->piDel = -1;

  getopts(cfg, argc, argv);
  srand(cfg->seed);
   
  if(cfg->probs || cfg->rates || cfg->hyper || cfg->parms || cfg->time || cfg->dic){
    process_mask = umask(0);
    mkdir(cfg->outDir, S_IRWXU | S_IRWXG | S_IRWXO); 
    umask(process_mask);
    chdir(cfg->outDir); 
  }
   
  if(cfg->probs){
    process_mask = umask(0);
    mkdir("probs", S_IRWXU | S_IRWXG | S_IRWXO);
    umask(process_mask);
  }
  
  if(cfg->rates){
    process_mask = umask(0);
    mkdir("rates", S_IRWXU | S_IRWXG | S_IRWXO);
    umask(process_mask);
  }
  
  if(cfg->hyper){
    process_mask = umask(0);
    mkdir("hyper", S_IRWXU | S_IRWXG | S_IRWXO);
    umask(process_mask);
  }
  
  if(cfg->parms){
    process_mask = umask(0);
    mkdir("parms", S_IRWXU | S_IRWXG | S_IRWXO); 
    umask(process_mask);
  }
  
  if(cfg->time){
    process_mask = umask(0);
    mkdir("time", S_IRWXU | S_IRWXG | S_IRWXO);
    umask(process_mask);
  }
  
  if(cfg->dic){
    process_mask = umask(0);
    mkdir("diagnostics", S_IRWXU | S_IRWXG | S_IRWXO);
    umask(process_mask);
  }
  
  return cfg;
}