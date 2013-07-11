#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void printRates(Chain *a, Config *cfg){

  int n, i, G = a->G, niter = cfg->M - cfg->burnin;
  num_t accD, accC, accPhi, accAlp, accDel, accEps;
  char file[BUF];
  FILE *fp;
  
  if(cfg->ratesFlag){
    fprintf(cfg->log, "  Printing Metropolis acceptance rates.\n");
  
    sprintf(file, "../out/rates/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "w"); 
    
    if(fp == NULL){
      printf("ERROR: unable to create file, %s\n", file);
      return;
    }
       
    fprintf(fp, "d c phi alp del meanEps\n");
    
    accD    = a->accD;
    accD   /= niter;
    
    accC    = a->accC[0];
    accC   /= niter;
    
    accPhi  = a->accPhi[0];
    accPhi /= niter;

    accAlp  = a->accAlp[0];
    accAlp /= niter;
    
    accDel  = a->accDel[0];
    accDel /= niter;
  
    accEps = 0;  
    for(n = 0; n < cfg->N; ++n)
      accEps += a->accEps[iG(n, 0)];
    accEps /= (niter * cfg->N);
    
    fprintf(fp, NUM_TF, accD);   fprintf(fp, " ");
    fprintf(fp, NUM_TF, accC);   fprintf(fp, " ");
    fprintf(fp, NUM_TF, accPhi); fprintf(fp, " ");
    fprintf(fp, NUM_TF, accAlp); fprintf(fp, " ");
    fprintf(fp, NUM_TF, accDel); fprintf(fp, " ");
    fprintf(fp, NUM_TF, accEps); fprintf(fp, " ");
    fprintf(fp, "\n");

    for(i = 1; i < cfg->N; ++i){
    
      accC    = a->accC[i];
      accC   /= niter;
    
      accPhi  = a->accPhi[i];
      accPhi /= niter;

      accAlp  = a->accAlp[i];
      accAlp /= niter;
    
      accDel  = a->accDel[i];
      accDel /= niter;
  
      accEps = 0;  
      for(n = 0; n < cfg->N; ++n)
        accEps += a->accEps[iG(n, i)];
      accEps /= (niter * cfg->N);
      
      fprintf(fp, ". ");
      fprintf(fp, NUM_TF, accC);   fprintf(fp, " ");
      fprintf(fp, NUM_TF, accPhi); fprintf(fp, " ");
      fprintf(fp, NUM_TF, accAlp); fprintf(fp, " ");
      fprintf(fp, NUM_TF, accDel); fprintf(fp, " ");
      fprintf(fp, NUM_TF, accEps); fprintf(fp, " ");
      fprintf(fp, "\n");
    }
    
    for(i = cfg->N; i < cfg->G; ++i){
    
      accPhi  = a->accPhi[i];
      accPhi /= niter;

      accAlp  = a->accAlp[i];
      accAlp /= niter;
    
      accDel  = a->accDel[i];
      accDel /= niter;
  
      accEps = 0;  
      
      for(n = 0; n < cfg->N; ++n)
        accEps += a->accEps[iG(n, i)];
      accEps /= (niter * cfg->N);
      
      fprintf(fp, ". . ");
      fprintf(fp, NUM_TF, accPhi); fprintf(fp, " ");
      fprintf(fp, NUM_TF, accAlp); fprintf(fp, " ");
      fprintf(fp, NUM_TF, accDel); fprintf(fp, " ");
      fprintf(fp, NUM_TF, accEps); fprintf(fp, " ");
      fprintf(fp, "\n");
    }
    
    fclose(fp);
  }
}