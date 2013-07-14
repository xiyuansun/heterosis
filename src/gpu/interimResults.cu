#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

__global__ void updateProbs(Chain *a){
  int g = IDX;

  if(a->m > a->burnin){
    if(g < a->G){
      a->dex[g] += ((a->alp[g] * a->alp[g]) > 1e-6);
  
      if(a->heterosis){
        a->hph[g] += (a->del[g] > fabs(a->alp[g]));
        a->lph[g] += (a->del[g] < -fabs(a->alp[g]));
        a->mph[g] += (fabs(a->del[g]) > 1e-6);
      }
    }
  }  
}

__global__ void updateM(Chain* a){
  ++a->m;
}

void sumLogLik_kernel(Chain *a){
  if(a->m > a->burnin)
    ++a->sumLogLik;
}

void interimResults(Chain *host_a, Chain *dev_a, Config *cfg){
  FILE *fp;
  char file[BUF];
  int n, g, G = cfg->G;
  num_t tmp, *tmpv;
  
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
      CUDA_CALL(cudaMemcpy(&(tmp), &(dev_a->sigC), sizeof(num_t), cudaMemcpyDeviceToHost));
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
  
  updateProbs<<<G_GRID, G_BLOCK>>>(dev_a);
  
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

    CUDA_CALL(cudaMemcpy(tmpv, host_a->phi, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));    
    for(g = 0; g < cfg->G; ++g){
      fprintf(fp, NUM_TF, tmpv[g]);
      fprintf(fp, " ");
    }
    
    CUDA_CALL(cudaMemcpy(tmpv, host_a->alp, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
    for(g = 0; g < cfg->G; ++g){
      fprintf(fp, NUM_TF, tmpv[g]);
      fprintf(fp, " ");
    }
    
    CUDA_CALL(cudaMemcpy(tmpv, host_a->del, cfg->N * sizeof(num_t), cudaMemcpyDeviceToHost));
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
  sumLogLik_kernel<<<1, 1>>>(dev_a);
  
  ++cfg->m;
  updateM<<<1, 1>>>(dev_a);
}