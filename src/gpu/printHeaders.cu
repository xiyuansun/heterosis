#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void printHeaders(Chain *host_a, Chain *dev_a, Config *cfg){
  FILE *fp;
  int n, g, G = cfg->G;
  char file[BUF];
  num_t tmp, *tmpv,;
 
  /* differential expression and heterosis probabilities */
  
  sprintf(file, "../out/probs/chain%d.txt", cfg->chainNum);
  fp = fopen(file, "w");
  
  if(fp == NULL){
    printf("ERROR: unable to create file, %s\n", file);
    return;
  }  
  
  fprintf(fp, "diff-expression high-parent-heterosis low-parent-heterosis mid-parent-heterosis\n");
  fclose(fp);
  
  /* hyperparameters */
  
  if(cfg->hyperFlag){
    sprintf(file, "../out/hyper/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "w");
    
    if(fp == NULL){
      printf("ERROR: unable to create file, %s\n", file);
      return;
    }
    
    fprintf(fp, "sigma-c d tau theta-phi theta-alpha theta-delta ");
    fprintf(fp, "sigma-phi sigma-alpha sigma-delta pi-alpha pi-delta\n");
    
    /* print initial values */
    
    CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->sigC), sizeof(num_t), cudaMemcpyDeviceToHost));
	fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
	
	CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->d), sizeof(num_t), cudaMemcpyDeviceToHost));
	fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
	
	CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->tau), sizeof(num_t), cudaMemcpyDeviceToHost));
	fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
	
	CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->thePhi), sizeof(num_t), cudaMemcpyDeviceToHost));
	fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
	
	CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->theAlp), sizeof(num_t), cudaMemcpyDeviceToHost));
	fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
	
	CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->theDel), sizeof(num_t), cudaMemcpyDeviceToHost));
	fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
	
	CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->sigPhi), sizeof(num_t), cudaMemcpyDeviceToHost));
	fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
	
	CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->sigAlp), sizeof(num_t), cudaMemcpyDeviceToHost));
	fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
	
	CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->sigDel), sizeof(num_t), cudaMemcpyDeviceToHost));
	fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
	
	CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->piAlp), sizeof(num_t), cudaMemcpyDeviceToHost));
	fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
	
	CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->piDel), sizeof(num_t), cudaMemcpyDeviceToHost));
	fprintf(fp, NUM_TF, tmp); fprintf(fp, "\n");
    
    fclose(fp);
  }
  
  /* parameters */
  
  if(cfg->parmsFlag){
    sprintf(file, "../out/parms/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "w");
    
    if(fp == NULL){
      printf("ERROR: unable to create file, %s\n", file);
      return;
    }
  
    /* print header */
    
    for(n = 0; n < cfg->N; ++n)
      fprintf(fp, "c%d ", n);    
     
    for(g = 0; g < cfg->G; ++g)
      fprintf(fp, "phi%d ", g);
    
    for(g = 0; g < cfg->G; ++g)
      fprintf(fp, "alpha%d ", g);
    
    for(g = 0; g < cfg->G; ++g)
      fprintf(fp, "delta%d ", g);
      
    for(g = 0; g < cfg->G; ++g)
      fprintf(fp, "eta%d ", g);
      
    for(g = 0; g < cfg->G; ++g)
      for(n = 0; n < cfg->N; ++n)
        fprintf(fp, "epsilon_lib%d_gene%d ", n, g);

    fprintf(fp, "\n");   
    
    /* print initial values */

    tmpv = (num_t*) malloc(cfg->N * cfg->G * sizeof(num_t));
  
    CUDA_CALL(cudaMemcpy(tmpv, host_a->c, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
    for(n = 0; n < cfg->N; ++n){
      fprintf(fp, NUM_TF, tmpv[n]);
      fprintf(fp, " ");
    }
    
    CUDA_CALL(cudaMemcpy(tmpv, host_a->phi, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
    for(g = 0; g < cfg->G; ++g){
      fprintf(fp, NUM_TF, tmpv[g]);
      fprintf(fp, " ");
    }
    
    CUDA_CALL(cudaMemcpy(tmpv, host_a->alp, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
    for(g = 0; g < cfg->G; ++g){
      fprintf(fp, NUM_TF, tmpv[g]);
      fprintf(fp, " ");
    }
    
    CUDA_CALL(cudaMemcpy(tmpv, host_a->del, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
    for(g = 0; g < cfg->G; ++g){
      fprintf(fp, NUM_TF, tmpv[g]);
      fprintf(fp, " ");
    }
    
    CUDA_CALL(cudaMemcpy(tmpv, host_a->eta, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
    for(g = 0; g < cfg->G; ++g){
      fprintf(fp, NUM_TF, tmpv[g]);
      fprintf(fp, " ");
    }
    
    CUDA_CALL(cudaMemcpy(tmpv, host_a->eps, cfg->N * cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
    for(g = 0; g < cfg->G; ++g)
      for(n = 0; n < cfg->N; ++n){
        fprintf(fp, NUM_TF, tmpv[iG(n, g)]);
        fprintf(fp, " ");
      }
     
    fprintf(fp, "\n");  
    fclose(fp);
    free(tmpv);
  }
  
  /* acceptance rates of Metropolis steps */
  
  if(cfg->ratesFlag){
  
    sprintf(file, "../out/rates/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "w"); 
    
    if(fp == NULL){
      printf("ERROR: unable to create file, %s\n", file);
      return;
    }
       
    fprintf(fp, "d c phi alp del meanEps\n");
    fclose(fp);
  }
  
  /* time spent sampling each parameter */
  
  if(cfg->timeFlag){
  
    sprintf(file, "../out/time/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "w"); 
  
    if(fp == NULL){
      printf("ERROR: unable to create file, %s\n", file);
      return;
    }
    
    fprintf(fp, "c tau pi-alpha pi-delta d theta-phi theta-alpha theta-delta sigma-c ");
    fprintf(fp, "sigma-phi sigma-alpha sigma-delta eta epsilon phi alpha delta\n");
    
    fclose(fp);  
  }
}