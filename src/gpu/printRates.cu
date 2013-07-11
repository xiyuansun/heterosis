#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void printRates(Chain *a, Config *cfg){

  int n, i, G = a->G, niter = cfg->M - cfg->burnin;
  num_t raccD, raccC, raccPhi, raccAlp, raccDel, raccEps;
  char file[BUF];
  FILE *fp;
  
  if(cfg->ratesFlag){
    sprintf(file, "../out/rates/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "w"); 
    
    if(fp == NULL){
      printf("ERROR: unable to create file, %s\n", file);
      return;
    }
       
    fprintf(fp, "d c phi alp del meanEps\n");
    
    raccD    = (num_t) a->accD;
    raccD   /= niter;
    
    raccC    = (num_t) a->accC[0];
    raccC   /= niter;
    
    raccPhi  = (num_t) a->accPhi[0];
    raccPhi /= niter;

    raccAlp  = (num_t) a->accAlp[0];
    raccAlp /= niter;
    
    raccDel  = (num_t) a->accDel[0];
    raccDel /= niter;
  
    raccEps = 0;  
    for(n = 0; n < cfg->N; ++n)
      raccEps += (num_t) a->accEps[iG(n, 0)];
    raccEps /= (niter * cfg->N);
    
    fprintf(fp, NUM_TF, raccD);   fprintf(fp, " ");
    fprintf(fp, NUM_TF, raccC);   fprintf(fp, " ");
    fprintf(fp, NUM_TF, raccPhi); fprintf(fp, " ");
    fprintf(fp, NUM_TF, raccAlp); fprintf(fp, " ");
    fprintf(fp, NUM_TF, raccDel); fprintf(fp, " ");
    fprintf(fp, NUM_TF, raccEps); fprintf(fp, " ");
    fprintf(fp, "\n");

    for(i = 1; i < cfg->N; ++i){
    
      raccC    = (num_t) a->accC[i];
      raccC   /= niter;
    
      raccPhi  = (num_t) a->accPhi[i];
      raccPhi /= niter;

      raccAlp  = (num_t) a->accAlp[i];
      raccAlp /= niter;
    
      raccDel  = (num_t) a->accDel[i];
      raccDel /= niter;
  
      raccEps = 0;  
      for(n = 0; n < cfg->N; ++n)
        raccEps += (num_t) a->accEps[iG(n, i)];
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
    
      raccPhi  = (num_t) a->accPhi[i];
      raccPhi /= niter;

      raccAlp  = (num_t) a->accAlp[i];
      raccAlp /= niter;
    
      raccDel  = (num_t) a->accDel[i];
      raccDel /= niter;
  
      raccEps = 0;  
      
      for(n = 0; n < cfg->N; ++n)
        raccEps += (num_t) a->accEps[iG(n, i)];
      raccEps /= (niter * cfg->N);
      
      fprintf(fp, ". . ");
      fprintf(fp, NUM_TF, raccPhi); fprintf(fp, " ");
      fprintf(fp, NUM_TF, raccAlp); fprintf(fp, " ");
      fprintf(fp, NUM_TF, raccDel); fprintf(fp, " ");
      fprintf(fp, NUM_TF, raccEps); fprintf(fp, " ");
      fprintf(fp, "\n");
    }
    
    fclose(fp);
  }
}