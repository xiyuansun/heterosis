#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void interimResults(Chain *host_a, Chain *dev_a, Config *cfg){
  FILE *fp;
  char file[BUF];
  int n, g, G = cfg->G;
  num_t tmp, *tmpv, *phi = NULL, *alp = NULL, *del = NULL;
  
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
      CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->sigC0), sizeof(num_t), cudaMemcpyDeviceToHost));
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    }
    
    if(cfg->constD){
      fprintf(fp, ". ");
    } else {
      CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->d), sizeof(num_t), cudaMemcpyDeviceToHost));
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    }
    
    if(cfg->constTau){
      fprintf(fp, ". ");
    } else {
      CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->tau), sizeof(num_t), cudaMemcpyDeviceToHost));
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    }
    
    if(cfg->constThePhi){
      fprintf(fp, ". ");
    } else {
      CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->thePhi), sizeof(num_t), cudaMemcpyDeviceToHost));
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    }
    
    if(cfg->constTheAlp){
      fprintf(fp, ". ");
    } else {
      CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->theAlp), sizeof(num_t), cudaMemcpyDeviceToHost));
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    }
    
    if(cfg->constTheDel){
      fprintf(fp, ". ");
    } else {
      CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->theDel), sizeof(num_t), cudaMemcpyDeviceToHost));
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    }
    
    if(cfg->constSigPhi){
      fprintf(fp, ". ");
    } else {
      CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->sigPhi), sizeof(num_t), cudaMemcpyDeviceToHost));
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    }
    
    if(cfg->constSigAlp){
      fprintf(fp, ". ");
    } else {
      CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->sigAlp), sizeof(num_t), cudaMemcpyDeviceToHost));
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    }
    
    if(cfg->constSigDel){
      fprintf(fp, ". ");
    } else {
      CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->sigDel), sizeof(num_t), cudaMemcpyDeviceToHost));
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    }
    
    if(cfg->constPiAlp){
      fprintf(fp, ". ");
    } else {
      CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->piAlp), sizeof(num_t), cudaMemcpyDeviceToHost));
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    }
    
    if(cfg->constPiDel){
      fprintf(fp, ". ");
    } else {
      CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->piDel), sizeof(num_t), cudaMemcpyDeviceToHost));
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    }

    fprintf(fp, "\n");
    
    fclose(fp);
  }
  
  /* heterosis and differential expression probabilities */
  
  phi  = (num_t*) malloc(cfg->G * sizeof(num_t));
  alp  = (num_t*) malloc(cfg->G * sizeof(num_t));
  del  = (num_t*) malloc(cfg->G * sizeof(num_t));  

  CUDA_CALL(cudaMemcpy(phi, host_a->phi, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(alp, host_a->alp, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(del, host_a->del, cfg->G * sizeof(num_t), cudaMemcpyDeviceToHost));
  
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
  
  /* parameters */
  
  if(cfg->parmsFlag){
    sprintf(file, "../out/parms/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "a");
    
    if(fp == NULL){
      printf("ERROR: unable to open file, %s\n", file);
      return;
    }
    
    tmpv = (num_t*) malloc(cfg->N * cfg->G * sizeof(num_t));
    
    CUDA_CALL(cudaMemcpy(tmpv, host_a->c, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
    for(n = 0; n < cfg->N; ++n){
      fprintf(fp, NUM_TF, tmpv[n]);
      fprintf(fp, " ");
    }
    
    for(g = 0; g < cfg->G; ++g){
      fprintf(fp, NUM_TF, phi[g]);
      fprintf(fp, " ");
    }
    
    for(g = 0; g < cfg->G; ++g){
      fprintf(fp, NUM_TF, alp[g]);
      fprintf(fp, " ");
    }
    
    for(g = 0; g < cfg->G; ++g){
      fprintf(fp, NUM_TF, del[g]);
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
  
  free(phi);
  free(alp);
  free(del);
  
  /* time spent in each sampler */
  
  if(cfg->timeFlag){
  
    sprintf(file, "../out/time/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "a"); 
  
    if(fp == NULL){
      printf("ERROR: unable to create file, %s\n", file);
      return;
    }
    
    CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->timeC), sizeof(num_t), cudaMemcpyDeviceToHost));
    fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");

    CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->timeTau), sizeof(num_t), cudaMemcpyDeviceToHost));
    fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");

    CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->timePiAlp), sizeof(num_t), cudaMemcpyDeviceToHost));
    fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");

    CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->timePiDel), sizeof(num_t), cudaMemcpyDeviceToHost));
    fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    
    CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->timeD), sizeof(num_t), cudaMemcpyDeviceToHost));
    fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    
    CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->timeThePhi), sizeof(num_t), cudaMemcpyDeviceToHost));
    fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    
    CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->timeTheAlp), sizeof(num_t), cudaMemcpyDeviceToHost));
    fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    
    CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->timeTheDel), sizeof(num_t), cudaMemcpyDeviceToHost));
    fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");

    CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->timeSigC), sizeof(num_t), cudaMemcpyDeviceToHost));
    fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");

    CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->timeSigPhi), sizeof(num_t), cudaMemcpyDeviceToHost));
    fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    
    CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->timeSigAlp), sizeof(num_t), cudaMemcpyDeviceToHost));
    fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    
    CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->timeSigDel), sizeof(num_t), cudaMemcpyDeviceToHost));
    fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
  
    CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->timeEta), sizeof(num_t), cudaMemcpyDeviceToHost));
    fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    
    CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->timeEps), sizeof(num_t), cudaMemcpyDeviceToHost));
    fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    
    CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->timephi), sizeof(num_t), cudaMemcpyDeviceToHost));
    fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    
    CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->timeAlp), sizeof(num_t), cudaMemcpyDeviceToHost));
    fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    
    CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->timeDel), sizeof(num_t), cudaMemcpyDeviceToHost));
    fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
    fprintf(fp, "\n");

    fclose(fp);  
  }
  
  ++cfg->m;
  ++a->m;
}