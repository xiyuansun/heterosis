#include <Chain.h>
#include <Config.h>
#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void printHyper(Chain *a, Config *cfg){

  int m;
  num_t tmp;
  num_t sigC, d, tau, thePhi, theAlp, theDel, sigPhi, sigAlp, sigDel, piAlp, piDel;
  char file[BUF];
  FILE *fp;
  
  if(cfg->hyperFlag){
    sprintf(file, "../out/hyper/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "w");
    
    if(fp == NULL){
      printf("ERROR: unable to create file, %s\n", file);
      return;
    }
    
    sigC   = (num_t*) malloc(cfg->M * sizeof(num_t));
    d      = (num_t*) malloc(cfg->M * sizeof(num_t));
    tau    = (num_t*) malloc(cfg->M * sizeof(num_t));
    thePhi = (num_t*) malloc(cfg->M * sizeof(num_t));
    theAlp = (num_t*) malloc(cfg->M * sizeof(num_t));
    theDel = (num_t*) malloc(cfg->M * sizeof(num_t));
    sigPhi = (num_t*) malloc(cfg->M * sizeof(num_t));
    sigAlp = (num_t*) malloc(cfg->M * sizeof(num_t));
    sigDel = (num_t*) malloc(cfg->M * sizeof(num_t));
    piAlp  = (num_t*) malloc(cfg->M * sizeof(num_t));
    piDel  = (num_t*) malloc(cfg->M * sizeof(num_t));
    
    CUDA_CALL(cudaMemcpy(sigC,   a->sigC,   cfg->M * sizeof(num_t), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(d,      a->d,      cfg->M * sizeof(num_t), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(tau,    a->tau,    cfg->M * sizeof(num_t), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(thePhi, a->thePhi, cfg->M * sizeof(num_t), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(theAlp, a->theAlp, cfg->M * sizeof(num_t), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(theDel, a->theDel, cfg->M * sizeof(num_t), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(sigPhi, a->sigPhi, cfg->M * sizeof(num_t), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(sigAlp, a->sigAlp, cfg->M * sizeof(num_t), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(sigDel, a->sigDel, cfg->M * sizeof(num_t), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(piAlp,  a->piAlp,  cfg->M * sizeof(num_t), cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(piDel,  a->piDel,  cfg->M * sizeof(num_t), cudaMemcpyDeviceToHost));
    
    fprintf(fp, "sigC d tau thePhi theAlp theDel sigPhi sigAlp sigDel piAlp piDel\n");
    
    for(m = 0; m <= cfg->M; ++m){ 
      fprintf(fp, NUM_TF, sigC[m]);   fprintf(fp, " ");   
      fprintf(fp, NUM_TF, d[m]);      fprintf(fp, " "); 
      fprintf(fp, NUM_TF, tau[m]);    fprintf(fp, " "); 
      fprintf(fp, NUM_TF, thePhi[m]); fprintf(fp, " ");  
      fprintf(fp, NUM_TF, theAlp[m]); fprintf(fp, " "); 
      fprintf(fp, NUM_TF, theDel[m]); fprintf(fp, " "); 
      fprintf(fp, NUM_TF, sigPhi[m]); fprintf(fp, " "); 
      fprintf(fp, NUM_TF, sigAlp[m]); fprintf(fp, " "); 
      fprintf(fp, NUM_TF, sigDel[m]); fprintf(fp, " "); 
      fprintf(fp, NUM_TF, piAlp[m]);  fprintf(fp, " "); 
      fprintf(fp, NUM_TF, piDel[m]);  fprintf(fp, " "); 
      
      fprintf(fp, "\n");
    }
    
    free(sigC);
    free(d);
    free(tau);
    free(thePhi);
    free(theAlp);
    free(theDel);
    free(sigPhi);
    free(sigAlp);
    free(sigDel);
    free(piAlp);
    free(piDel);
    
    fclose(fp);
  }
}