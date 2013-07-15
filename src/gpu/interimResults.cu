#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

__global__ void updateProbs(Chain *a, int heterosis){
  int g = IDX;

  if(a->m > a->burnin){
    if(g < a->G){
      a->dex[g] += ((a->alp[g] * a->alp[g]) > 1e-6);
  
      if(heterosis){
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

__host__ void printHyper(Chain *host_a, Chain *dev_a, Config *cfg){
  
  char file[BUF];
  FILE *fp;
  num_t tmp;
  
  if(cfg->hyper){
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
}


__host__ void printParms(Chain *host_a, Chain *dev_a, Config *cfg){

  int n, g, G = cfg->G;
  char file[BUF];
  FILE *fp;
  num_t tmp, *tmpv;

  if(cfg->parms){
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
    for(n = 0; n < cfg->N; ++n)
      for(g = 0; g < cfg->G; ++g){
        fprintf(fp, NUM_TF, tmpv[iG(n, g)]);
        fprintf(fp, " ");
      }
      
    fprintf(fp, "\n"); 
        
    fclose(fp);
    free(tmpv);
  }
} 

__host__ void printTime(Chain *host_a, Chain *dev_a, Config *cfg){
  char file[BUF];
  FILE *fp;

  if(cfg->time){
  
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
}

__host__ void interimResults(Chain *host_a, Chain *dev_a, Config *cfg){

  if(cfg->probs) 
    updateProbs<<<G_GRID, G_BLOCK>>>(dev_a, cfg->heterosis);
  
  printHyper(host_a, dev_a, cfg);
  printParms(host_a, dev_a, cfg);
  printTime(host_a, dev_a, cfg);
  
  if(cfg->dic)
    updateDICprep(dev_a, cfg);
  
  ++cfg->m;
  updateM<<<1, 1>>>(dev_a);
}
