#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

void summarizeChain(Chain *a, Config *cfg){
  int n, g, i, N = cfg->N, G = cfg->G,  niter = cfg->M - cfg->burnin;
  num_t accD, *accC, *accPhi, *accAlp, *accDel, *accEps; 
  FILE *fp;
  char file[BUF];

  /* differential expression and heterosis probabilities */
  
  if(cfg->verbose){
    printf("  Printing differential expression ");
    
    if(cfg->heterosis)
      printf("and heterosis ");
      
    printf("probabilities.\n");
  }
  
  sprintf(file, "../out/probs/chain%d.txt", cfg->chainNum);
  fp = fopen(file, "a");
  
  if(fp == NULL){
    printf("ERROR: unable to create file, %s\n", file);
    return;
  }
  
  for(g = 0; g < G; ++g){
    fprintf(fp, NUM_TF, ((num_t) a->dex[g]) / niter); fprintf(fp, " ");
    
    if(cfg->heterosis){
      fprintf(fp, NUM_TF, ((num_t) a->hph[g]) / niter); fprintf(fp, " ");
      fprintf(fp, NUM_TF, ((num_t) a->lph[g]) / niter); fprintf(fp, " ");
      fprintf(fp, NUM_TF, ((num_t) a->mph[g]) / niter); fprintf(fp, " ");
    }
    fprintf(fp, "\n");
  }
  
  fclose(fp);
  
  if(cfg->ratesFlag){
  
    if(cfg->verbose)
      printf("  Printing acceptance rates of metropolis steps.\n");
  
    sprintf(file, "../out/rates/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "a"); 
    
    if(fp == NULL){
      printf("ERROR: unable to create file, %s\n", file);
      return;
    }
        
    accC = (num_t*) malloc(cfg->N * sizeof(num_t));
    accPhi = (num_t*) malloc(cfg->G * sizeof(num_t));
    accAlp = (num_t*) malloc(cfg->G * sizeof(num_t));
    accDel = (num_t*) malloc(cfg->G * sizeof(num_t));
    accEps = (num_t*) malloc(cfg->G * sizeof(num_t));   
      
    accD = a->accD / (num_t) niter;
    
    for(n = 0; n < cfg->N; ++n)
      accC[n] = a->accC[n] / (num_t) niter;
      
    for(g = 0; g < cfg->G; ++g){
      accPhi[g] = a->accPhi[g] / (num_t) niter;
      accAlp[g] = a->accAlp[g] / (num_t) niter;
      accDel[g] = a->accDel[g] / (num_t) niter;
      
      accEps[g] = 0;
      for(n = 0; n < cfg->N; ++n)
        accEps[g] += a->accEps[iG(n, g)] / (num_t) niter;
    }  
    
    for(i = 0; i < MAX_NG; ++i){
      if(!i){
        fprintf(fp, NUM_TF, accD);   fprintf(fp, " ");
      } else {
        fprintf(fp, ". ");
      }
    
      if(i < cfg->N){
        fprintf(fp, NUM_TF, accC[g]);   fprintf(fp, " ");
      } else {
        fprintf(fp, ". ");
      }
      
      if(i < cfg->G){
        fprintf(fp, NUM_TF, accPhi[g]); fprintf(fp, " ");
        fprintf(fp, NUM_TF, accAlp[g]); fprintf(fp, " ");
        fprintf(fp, NUM_TF, accDel[g]); fprintf(fp, " ");
        fprintf(fp, NUM_TF, accEps[g]); fprintf(fp, " ");
      } else {
        fprintf(fp, ". . . . ");
      }

      fprintf(fp, "\n");
    }
    
    free(accC);
    free(accPhi);
    free(accAlp);
    free(accDel);
    free(accEps);
    
    fclose(fp);
  }
}