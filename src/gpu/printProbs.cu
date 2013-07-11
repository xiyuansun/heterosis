#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

__host__ void printProbs(Chain *a, Config *cfg){
  int m, g, G = a->G, niter = cfg->M - cfg->burnin;
  num_t alp, del;
  num_t prob_de, prob_hph, prob_lph, prob_mph;
  FILE *fp;
  
  if(cfg->probsFlag){
    fp = fopen(cfg->probsFile, "w"); 
    
    if(fp == NULL){
      printf("ERROR: unable to create file, %s\n", cfg->probsFile);
      return;
    }
    
    fprintf(fp, "de ");
    if(cfg->heterosis)
      fprintf(fp, "hph lph mph");
    fprintf(fp, "\n");
    
    for(g = 0; g < cfg->G; ++g){
      prob_de = 0;
      
      for(m = cfg->burnin + 1; m <= cfg->M; ++m){
        alp = a->alp[iG(m, g)];
        prob_de += ((alp * alp) > 1e-6);
      }
      
      prob_de /= niter;
      fprintf(fp, NUM_TF, prob_de);
      fprintf(fp, " ");
      
      if(cfg->heterosis){ 
      
        prob_hph = 0;
        prob_lph = 0;
        prob_mph = 0;
      
        for(m = cfg->burnin + 1; m <= cfg->M; ++m){
          alp = a->alp[iG(m, g)];
          del = a->del[iG(m, g)];
          
          prob_hph += (del > fabs(alp));
          prob_lph += (del < -fabs(alp));
          prob_mph += (fabs(del) > 1e-6);
        }
      
        prob_hph /= niter;
        prob_lph /= niter;
        prob_mph /= niter;
        
        fprintf(fp, NUM_TF, prob_hph); fprintf(fp, " ");
        fprintf(fp, NUM_TF, prob_lph); fprintf(fp, " ");
        fprintf(fp, NUM_TF, prob_mph); fprintf(fp, " ");
      }
      fprintf(fp, "\n");
    }
    
    fclose(fp);
  }
}