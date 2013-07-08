#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

void printRates(Chain *a, Config *cfg){

  int n, i, niter = cfg->M - cfg->burnin;
  int M = cfg->M, N = cfg->N, G = cfg->G;
  num_t accD, accC, accPhi, accAlp, accDel, accEps;
  FILE *fp;
  
  if(cfg->ratesFlag){
    fp = fopen(cfg->ratesFile, "w");  
    
    if(fp == NULL){
      printf("ERROR: unable to create file, %s\n", cfg->ratesFile);
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
      accEps += a->accEps[iNG(n, 0)];
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
        accEps += a->accEps[iNG(n, i)];
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
        accEps += a->accEps[iNG(n, i)];
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