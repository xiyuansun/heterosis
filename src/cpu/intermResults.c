#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void intermResults(Chain *a, Config *cfg){
  FILE *fp;
  int n, g, N = cfg->N, G = cfg->G;
  
  /* hyperparameters */
  
  if(cfg->hyperFlag){
    sprintf(file, "../out/hyper/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "w+");
    
    if(fp == NULL){
      printf("ERROR: unable to open file, %s\n", file);
      return;
    }
    
    if(cfg->constSigC){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->sigC); fprintf(" ");
    }
    
    if(cfg->constD){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->d); fprintf(" ");
    }
    
    if(cfg->constTau){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->tau); fprintf(" ");
    }
    
    if(cfg->constThePhi){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->thePhi); fprintf(" ");
    }
    
    if(cfg->constTheAlp){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->theAlp); fprintf(" ");
    }
    
    if(cfg->constTheDel){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->theDel); fprintf(" ");
    }
    
    if(cfg->constSigPhi){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->sigPhi); fprintf(" ");
    }
    
    if(cfg->constSigAlp){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->sigAlp); fprintf(" ");
    }
    
    if(cfg->constSigDel){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->sigDel); fprintf(" ");
    }
    
    if(cfg->constPiAlp){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->sigPiAlp); fprintf(" ");
    }
    
    if(cfg->constPiDel){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->sigPiDel); fprintf(" ");
    }

    fprintf(fp, "\n");
    
    fclose(fp);
  }
  
  /* parameters */
  
  if(cfg->parmsFlag){
    sprintf(file, "../out/parms/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "w+");
    
    if(fp == NULL){
      printf("ERROR: unable to open file, %s\n", file);
      return;
    }
    
    for(n = 0; n < cfg->N; ++n){
      fprintf(fp, NUM_TF, a->c[n]);
      fprintf(fp, " ");
    }
    
    for(g = 0; g < cfg->G; ++g){
      fprintf(fp, NUM_TF, a->phi[g]);
      fprintf(fp, " ");
    }
    
    for(g = 0; g < cfg->G; ++g){
      fprintf(fp, NUM_TF, a->alp[g]);
      fprintf(fp, " ");
    }
    
    for(g = 0; g < cfg->G; ++g){
      fprintf(fp, NUM_TF, a->del[g]);
      fprintf(fp, " ");
    }
    
    for(g = 0; g < cfg->G; ++g){
      fprintf(fp, NUM_TF, a->eta[g]);
      fprintf(fp, " ");
    }
    
    for(g = 0; g < cfg->G; ++g)
      for(n = 0; n < cfg->N; ++n){
        fprintf(fp, NUM_TF, a->eps[iG(n, g)]);
        fprintf(fp, " ");
      }
      
    fprintf(fp, "\n");    
    fclose(fp);
  }
  
  if(a->m > cfg->burnin){
    for(g = 0; g < G; ++g){
      a->dex += ((alp[g] * alp[g]) > 1e-6);
  
      if(cfg->heterosis){
        a->hph += (del[g] > fabs(alp[g]));
        a->lph += (del[g] < -fabs(alp[g]));
        a->mph += (fabs(del[g]) > 1e-6);
      }
    }
  }
  
  /* time spent in each sampler */
  
  if(cfg->timeFlag){
  
    sprintf(file, "../out/time/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "w"); 
  
    if(fp == NULL){
      printf("ERROR: unable to create file, %s\n", file);
      return;
    }
    
    fprintf(fp, NUM_TF, cfg->timeC); fprintf(fp, " ");
    fprintf(fp, NUM_TF, cfg->timeTau); fprintf(fp, " ");
    fprintf(fp, NUM_TF, cfg->timePiAlp); fprintf(fp, " ");
    fprintf(fp, NUM_TF, cfg->timePiDel); fprintf(fp, " ");
    fprintf(fp, NUM_TF, cfg->timeD); fprintf(fp, " ");
    fprintf(fp, NUM_TF, cfg->timeThePhi); fprintf(fp, " ");
    fprintf(fp, NUM_TF, cfg->timeTheAlp); fprintf(fp, " ");
    fprintf(fp, NUM_TF, cfg->timeTheDel); fprintf(fp, " ");
    fprintf(fp, NUM_TF, cfg->timeSigC); fprintf(fp, " ");
    fprintf(fp, NUM_TF, cfg->timeSigPhi); fprintf(fp, " ");
    fprintf(fp, NUM_TF, cfg->timeSigAlp); fprintf(fp, " ");
    fprintf(fp, NUM_TF, cfg->timeSigDel); fprintf(fp, " ");
    fprintf(fp, NUM_TF, cfg->timeEta); fprintf(fp, " ");
    fprintf(fp, NUM_TF, cfg->timeEps); fprintf(fp, " ");
    fprintf(fp, NUM_TF, cfg->timePhi); fprintf(fp, " ");
    fprintf(fp, NUM_TF, cfg->timeAlp); fprintf(fp, " ");
    fprintf(fp, NUM_TF, cfg->timeDel); fprintf(fp, " ");
    fprintf(fp, "\n");

    fclose(fp);  
  }
  
  ++cfg->m;
  ++a->m;
}