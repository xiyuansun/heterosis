#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

void summarizeChain(Chain *a, Config *cfg){
  int n, g, i, N = cfg->N, G = cfg->G,  niter = cfg->M - cfg->burnin;
  num_t accD, accC, accPhi, accAlp, accDel, accEps; 

  /* differential expression and heterosis probabilities */
  
  sprintf(file, "../out/probs/chain%d.txt", cfg->chainNum);
  fp = fopen(file, "w");
  
  for(g = 0; g < G; ++g){
    
    fprintf(fp, NUM_TF, ((num_t) a->dex) / niter); printf(" ");
    
    if(cfg->heterosis){
      fprintf(fp, NUM_TF, ((num_t) a->hph) / niter); printf(" ");
      fprintf(fp, NUM_TF, ((num_t) a->lph) / niter); printf(" ");
      fprintf(fp, NUM_TF, ((num_t) a->mph) / niter); printf(" ");
    }
  }
 
  fclose(fp);
  
  /* acceptance rates of Metropolis steps */
  
  if(cfg->ratesFlag){
  
    sprintf(file, "../out/rates/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "w"); 
    
    if(fp == NULL){
      printf("ERROR: unable to create file, %s\n", file);
      return;
    }
       
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