#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void printProbs(Chain *a, Config *cfg){
  int m, g, G = cfg->G, niter = cfg->M - cfg->burnin;
  num_t *alp, *del;
  num_t prob_de, prob_hph, prob_lph, prob_mph;
  char file[BUF] = "../out/probs/chain";
  FILE *fp;
  
  sprintf(file, "../out/probs/chain%d.txt", cfg->chainNum);
  fp = fopen(file, "w"); 
    
  if(fp == NULL){
    printf("ERROR: unable to create file, %s\n", file);
    return;
  }
    
  fprintf(fp, "de ");
  if(cfg->heterosis)
    fprintf(fp, "hph lph mph");
  fprintf(fp, "\n");
    
  alp = (num_t*) malloc(cfg->N * cfg->G * sizeof(num_t));    
  del = (num_t*) malloc(cfg->N * cfg->G * sizeof(num_t));

  CUDA_CALL(cudaMemcpy(alp, a->alp, cfg->N * cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost);
  CUDA_CALL(cudaMemcpy(del, a->del, cfg->N * cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost);
    
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
        prob_hph += (del[iG(m, g)] > fabs(alp[iG(m, g)]));
        prob_lph += (del[iG(m, g)] < -fabs(alp[iG(m, g)]));
        prob_mph += (fabs(del[iG(m, g)]) > 1e-6);
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

  free(alp);
  free(del);   
  fclose(fp);
}