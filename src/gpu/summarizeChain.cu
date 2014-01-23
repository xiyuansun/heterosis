#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

__host__ void printProbs(Chain *host_a, Chain *dev_a, Config *cfg){

  int g, G = cfg->G,  niter = cfg->M - cfg->burnin, *dex, *hph, *lph, *mph;
  FILE *fp;
  char file[BUF];

  if(cfg->probs){
    if(cfg->verbose){
	  printf("  Printing differential expression ");
	
	  if(cfg->heterosis)
		printf("and heterosis ");
	  
	  printf("probabilities.\n");
	}
  
	sprintf(file, "probs-chain%d.txt", cfg->chainNum);
	fp = fopen(file, "a");
  
	if(fp == NULL){
	  fprintf(stderr, "ERROR: unable to create file, %s\n", file);
	  exit(EXIT_FAILURE);
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
  }
}

__host__ void printRates(Chain *host_a, Chain *dev_a, Config *cfg){

  int n, g, i, N = cfg->N, G = cfg->G,  niter = cfg->M - cfg->burnin;
  int nAccD, *nAccC, *nAccPhi, *nAccAlp, *nAccDel, *nAccEps;
  num_t accD, *accC, *accPhi, *accAlp, *accDel, *accEps;
  FILE *fp;
  char file[BUF];

  if(cfg->rates){
  
    if(cfg->verbose)
      printf("  Printing acceptance rates of metropolis steps.\n");
  
    sprintf(file, "rates-chain%d.txt", cfg->chainNum);
    fp = fopen(file, "a"); 
    
    if(fp == NULL){
      fprintf(stderr, "ERROR: unable to create file, %s\n", file);
      exit(EXIT_FAILURE);
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
        accEps[g] += nAccEps[iG(n, g)];
     
       accEps[g] /= (num_t) (niter * cfg->N);
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

__host__ void printDIC(Chain *host_a, Chain *dev_a, Config *cfg){

  FILE *fp;
  num_t tmp;
  
  if(cfg->dic){
  
	fp = fopen("dic.txt", "a");
  
	if(fp == NULL){
	  fprintf(stderr, "ERROR: unable to create file, dic.txt.\n");
	  exit(EXIT_FAILURE);
	}
  
	dic<<<1, 1>>>(dev_a);
	CUDA_CALL(cudaMemcpy(&tmp, &(dev_a->dic), sizeof(num_t), cudaMemcpyDeviceToHost));
	
	fprintf(fp, NUM_TF, tmp);
	fprintf(fp, "\n");
	fclose(fp);
  }  
}

__host__ void summarizeChain(Chain *host_a, Chain *dev_a, Config *cfg){

  printProbs(host_a, dev_a, cfg);
  printRates(host_a, dev_a, cfg);
  printDIC(host_a, dev_a, cfg);
}
