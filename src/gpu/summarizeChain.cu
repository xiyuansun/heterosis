#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

void summarizeChain(Chain *host_a, Chain *dev_a, Config *cfg){
  int n, g, i, G = cfg->G,  niter = cfg->M - cfg->burnin, *dex, *hph, *lph, *mph;
  num_t accD, *accC, *accPhi, *accAlp, *accDel, *accEps; 
  int nAccD, *nAccC, *nAccPhi, *nAccAlp, *nAccDel, *nAccEps;
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
  
  dex = (int*) malloc(cfg->G * sizeof(int));
  hph = (int*) malloc(cfg->G * sizeof(int));
  lph = (int*) malloc(cfg->G * sizeof(int));
  mph = (int*) malloc(cfg->G * sizeof(int));
  
  CUDA_CALL(cudaMemcpy(dex, host_a->dex, cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(hph, host_a->hph, cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(lph, host_a->lph, cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(mph, host_a->mph, cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
  
  for(g = 0; g < G; ++g){
    fprintf(fp, NUM_TF, ((num_t) dex[g]) / niter); fprintf(fp, " ");
    
    if(cfg->heterosis){
      fprintf(fp, NUM_TF, ((num_t) hph[g]) / niter); fprintf(fp, " ");
      fprintf(fp, NUM_TF, ((num_t) lph[g]) / niter); fprintf(fp, " ");
      fprintf(fp, NUM_TF, ((num_t) mph[g]) / niter); fprintf(fp, " ");
    }
    fprintf(fp, "\n");
  }
  
  free(dex);
  free(hph);
  free(lph);
  free(mph);
  fclose(fp);
  
  /* acceptance rates of Metropolis steps */
  
  if(cfg->ratesFlag){
  
    if(cfg->verbose)
      printf("  Printing acceptance rates of metropolis steps.\n");
  
    sprintf(file, "../out/rates/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "a"); 
    
    if(fp == NULL){
      printf("ERROR: unable to create file, %s\n", file);
      return;
    }
       
    nAccC = (int*) malloc(cfg->N * sizeof(int));
    nAccPhi = (int*) malloc(cfg->G * sizeof(int));
    nAccAlp = (int*) malloc(cfg->G * sizeof(int));
    nAccDel = (int*) malloc(cfg->G * sizeof(int));
    nAccEps = (int*) malloc(cfg->N * cfg->G * sizeof(int));   
       
    CUDA_CALL(cudaMemcpy(&(nAccD), &(dev_a->accD), sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(nAccC, host_a->accC, cfg->N * sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(nAccPhi, host_a->accPhi, cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(nAccAlp, host_a->accAlp, cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(nAccDel, host_a->accDel, cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(nAccEps, host_a->accEps, cfg->N * cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
        
    accC = (num_t*) malloc(cfg->N * sizeof(num_t));
    accPhi = (num_t*) malloc(cfg->G * sizeof(num_t));
    accAlp = (num_t*) malloc(cfg->G * sizeof(num_t));
    accDel = (num_t*) malloc(cfg->G * sizeof(num_t));
    accEps = (num_t*) malloc(cfg->G * sizeof(num_t));   
      
    accD = nAccD / (num_t) niter;
    
    for(n = 0; n < cfg->N; ++n)
      accC[n] = nAccC[n] / (num_t) niter;
      
    for(g = 0; g < cfg->G; ++g){
      accPhi[g] = nAccPhi[g] / (num_t) niter;
      accAlp[g] = nAccAlp[g] / (num_t) niter;
      accDel[g] = nAccDel[g] / (num_t) niter;
      
      accEps[g] = 0;
      for(n = 0; n < cfg->N; ++n)
        accEps[g] += nAccEps[iG(n, g)] / (num_t) niter;
    }  
    
    for(i = 0; i < MAX_NG; ++i){
      if(!i){
        fprintf(fp, NUM_TF, accD);   fprintf(fp, " ");
      } else {
        fprintf(fp, ". ");
      }
    
      if(i < cfg->N){
        fprintf(fp, NUM_TF, accC[i]);   fprintf(fp, " ");
      } else {
        fprintf(fp, ". ");
      }
      
      if(i < cfg->G){
        fprintf(fp, NUM_TF, accPhi[i]); fprintf(fp, " ");
        fprintf(fp, NUM_TF, accAlp[i]); fprintf(fp, " ");
        fprintf(fp, NUM_TF, accDel[i]); fprintf(fp, " ");
        fprintf(fp, NUM_TF, accEps[i]); fprintf(fp, " ");
      } else {
        fprintf(fp, ". . . . "); 
      }

      fprintf(fp, "\n");
    }
    
    free(nAccC);
    free(nAccPhi);
    free(nAccAlp);
    free(nAccDel);
    free(nAccEps);
    
    free(accC);
    free(accPhi);
    free(accAlp);
    free(accDel);
    free(accEps);
    
    fclose(fp);
  }
}