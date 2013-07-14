#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void sumLogLik_kernel1(Chain *a){ /* kernel<<<G, 1>>> */
  int n, g, N = a->N, G = a->G;
  num_t lam;
  
  for(g = 0; g < G; ++g){
    a->tmp1[g] = 0;
    for(n = 0; n < N; ++n){
      lam = a->c[n] + a->eps[iG(n, g)] 
          + mu(a, n, a->phi[g], a->alp[g], a->del[g]);
      a->tmp1[g] += a->y[iG(n, g)] * lam - exp(lam);
    }
  }
}

void sumLogLik_kernel2(Chain *a){
  a->sumLogLik += a->s1;
}

void sumLogLik(Chain *a, Config *cfg){
  int g;
  
  if(cfg->m <= cfg->burnin)
    return;

  sumLogLik_kernel1(a);
  
  /* pairwise sum in thrust */
  a->s1 = 0;
  for(g = 0; g < cfg->G; ++g)
    a->s1 += a->tmp1[g];
  
  sumLogLik_kernel2(a);
}

void intermResults(Chain *a, Config *cfg){
  FILE *fp;
  char file[BUF];
  int n, g, G = cfg->G;
  
  /* hyperparameters */
  
  if(cfg->hyperFlag){
    sprintf(file, "../out/hyper/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "a");
    
    if(fp == NULL){
      printf("ERROR: unable to open file, %s\n", file);
      return;
    }
    
    if(cfg->constSigC){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->sigC); fprintf(fp, " ");
    }
    
    if(cfg->constD){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->d); fprintf(fp, " ");
    }
    
    if(cfg->constTau){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->tau); fprintf(fp, " ");
    }
    
    if(cfg->constThePhi){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->thePhi); fprintf(fp, " ");
    }
    
    if(cfg->constTheAlp){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->theAlp); fprintf(fp, " ");
    }
    
    if(cfg->constTheDel){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->theDel); fprintf(fp, " ");
    }
    
    if(cfg->constSigPhi){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->sigPhi); fprintf(fp, " ");
    }
    
    if(cfg->constSigAlp){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->sigAlp); fprintf(fp, " ");
    }
    
    if(cfg->constSigDel){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->sigDel); fprintf(fp, " ");
    }
    
    if(cfg->constPiAlp){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->piAlp); fprintf(fp, " ");
    }
    
    if(cfg->constPiDel){
      fprintf(fp, ". ");
    } else {
      fprintf(fp, NUM_TF, a->piDel); fprintf(fp, " ");
    }

    fprintf(fp, "\n");
    
    fclose(fp);
  }
  
  /* parameters */
  
  if(cfg->parmsFlag){
    sprintf(file, "../out/parms/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "a");
    
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
      a->dex[g] += ((a->alp[g] * a->alp[g]) > 1e-6);
  
      if(cfg->heterosis){
        a->hph[g] += (a->del[g] > fabs(a->alp[g]));
        a->lph[g] += (a->del[g] < -fabs(a->alp[g]));
        a->mph[g] += (fabs(a->del[g]) > 1e-6);
      }
    }
  }
  
  /* time spent in each sampler */
  
  if(cfg->timeFlag){
  
    sprintf(file, "../out/time/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "a"); 
  
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
  
  /* update across-chain sum of model likelihoods */
  sumLogLik(a, cfg);
  
  ++cfg->m;
  ++a->m;
}