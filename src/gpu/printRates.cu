#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void printRates(Chain *host_a, Chain *dev_a, Config *cfg){

  int n, i, G = cfg->G, niter = cfg->M - cfg->burnin;
  int accD, *accC, *accPhi, *accAlp, *accDel, *accEps;
  num_t raccD, raccC, raccPhi, raccAlp, raccDel, raccEps;
  char file[BUF];
  FILE *fp;
    
  if(cfg->ratesFlag){
  
    float myTime;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
  
    fprintf(cfg->log, "  Printing Metropolis acceptance rates.\n");
  
    sprintf(file, "../out/rates/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "w"); 
    
    if(fp == NULL){
      printf("ERROR: unable to create file, %s\n", file);
      return;
    }
       
    fprintf(fp, "d c phi alp del meanEps\n");
    
    accC = (int*) malloc(cfg->N * sizeof(int));
    accPhi = (int*) malloc(cfg->G * sizeof(int));
    accAlp = (int*) malloc(cfg->G * sizeof(int));
    accDel = (int*) malloc(cfg->G * sizeof(int));
    accEps = (int*) malloc(cfg->N * cfg->G * sizeof(int));
    
    CUDA_CALL(cudaMemcpy(&(accD), &(dev_a->accD), sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(accC, host_a->accC, cfg->N * sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(accPhi, host_a->accPhi, cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(accAlp, host_a->accAlp, cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(accDel, host_a->accDel, cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(accEps, host_a->accEps, cfg->N * cfg->G * sizeof(int), cudaMemcpyDeviceToHost));
    
    raccD    = (num_t) accD;
    raccD   /= niter;
    
    raccC    = (num_t) accC[0];
    raccC   /= niter;
    
    raccPhi  = (num_t) accPhi[0];
    raccPhi /= niter;

    raccAlp  = (num_t) accAlp[0];
    raccAlp /= niter;
    
    raccDel  = (num_t) accDel[0];
    raccDel /= niter;
  
    raccEps = 0;  
    for(n = 0; n < cfg->N; ++n)
      raccEps += (num_t) accEps[iG(n, 0)];
    raccEps /= (niter * cfg->N);
    
    fprintf(fp, NUM_TF, raccD);   fprintf(fp, " ");
    fprintf(fp, NUM_TF, raccC);   fprintf(fp, " ");
    fprintf(fp, NUM_TF, raccPhi); fprintf(fp, " ");
    fprintf(fp, NUM_TF, raccAlp); fprintf(fp, " ");
    fprintf(fp, NUM_TF, raccDel); fprintf(fp, " ");
    fprintf(fp, NUM_TF, raccEps); fprintf(fp, " ");
    fprintf(fp, "\n");

    for(i = 1; i < cfg->N; ++i){
    
      raccC    = (num_t) accC[i];
      raccC   /= niter;
    
      raccPhi  = (num_t) accPhi[i];
      raccPhi /= niter;

      raccAlp  = (num_t) accAlp[i];
      raccAlp /= niter;
    
      raccDel  = (num_t) accDel[i];
      raccDel /= niter;
  
      raccEps = 0;  
      for(n = 0; n < cfg->N; ++n)
        raccEps += (num_t) accEps[iG(n, i)];
      raccEps /= (niter * cfg->N);
      
      fprintf(fp, ". ");
      fprintf(fp, NUM_TF, raccC);   fprintf(fp, " ");
      fprintf(fp, NUM_TF, raccPhi); fprintf(fp, " ");
      fprintf(fp, NUM_TF, raccAlp); fprintf(fp, " ");
      fprintf(fp, NUM_TF, raccDel); fprintf(fp, " ");
      fprintf(fp, NUM_TF, raccEps); fprintf(fp, " ");
      fprintf(fp, "\n");
    }
    
    for(i = cfg->N; i < cfg->G; ++i){
    
      raccPhi  = (num_t) accPhi[i];
      raccPhi /= niter;

      raccAlp  = (num_t) accAlp[i];
      raccAlp /= niter;
    
      raccDel  = (num_t) accDel[i];
      raccDel /= niter;
  
      raccEps = 0;  
      
      for(n = 0; n < cfg->N; ++n)
        raccEps += (num_t) accEps[iG(n, i)];
      raccEps /= (niter * cfg->N);
      
      fprintf(fp, ". . ");
      fprintf(fp, NUM_TF, raccPhi); fprintf(fp, " ");
      fprintf(fp, NUM_TF, raccAlp); fprintf(fp, " ");
      fprintf(fp, NUM_TF, raccDel); fprintf(fp, " ");
      fprintf(fp, NUM_TF, raccEps); fprintf(fp, " ");
      fprintf(fp, "\n");
    }

    free(accC);
    free(accPhi);
    free(accAlp);
    free(accDel);
    free(accEps);
    fclose(fp);
    
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&myTime, start, stop);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    
    fprintf(cfg->timeConfig, "%0.3f ", myTime/MILLISECS); /* elapsed time */
  }
}