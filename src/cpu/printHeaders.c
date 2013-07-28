#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/stat.h>
#include <unistd.h>


void printHeaders(Chain *a, Config *cfg){
  FILE *fp;
  int n, g, G = cfg->G;
  char file[BUF];
   
  chdir(cfg->outDir); 
   
  /* differential expression and heterosis probabilities */
  
  if(cfg->probs){
	sprintf(file, "probs-chain%d.txt", cfg->chainNum);
	fp = fopen(file, "w");
  
	if(fp == NULL){
	  fprintf(stderr, "ERROR: unable to create file, %s\n", file);
	  exit(EXIT_FAILURE);
	}  
  
	fprintf(fp, "diff-expression high-parent-heterosis low-parent-heterosis mid-parent-heterosis\n");
	fclose(fp);
  }
  
  /* hyperparameters */
  
  if(cfg->hyper){
    sprintf(file, "hyper-chain%d.txt", cfg->chainNum);
    fp = fopen(file, "w");
    
    if(fp == NULL){
      fprintf(stderr, "ERROR: unable to create file, %s\n", file);
      exit(EXIT_FAILURE);
    }
    
    fprintf(fp, "sigma-c d tau theta-phi theta-alpha theta-delta ");
    fprintf(fp, "sigma-phi sigma-alpha sigma-delta pi-alpha pi-delta\n");
    
    /* print initial values */
    
	fprintf(fp, NUM_TF, a->sigC); fprintf(fp, " ");
	fprintf(fp, NUM_TF, a->d); fprintf(fp, " ");
	fprintf(fp, NUM_TF, a->tau); fprintf(fp, " ");
	fprintf(fp, NUM_TF, a->thePhi); fprintf(fp, " ");
	fprintf(fp, NUM_TF, a->theAlp); fprintf(fp, " ");
	fprintf(fp, NUM_TF, a->theDel); fprintf(fp, " ");
	fprintf(fp, NUM_TF, a->sigPhi); fprintf(fp, " ");
	fprintf(fp, NUM_TF, a->sigAlp); fprintf(fp, " ");
	fprintf(fp, NUM_TF, a->sigDel); fprintf(fp, " ");
	fprintf(fp, NUM_TF, a->piAlp); fprintf(fp, " ");
	fprintf(fp, NUM_TF, a->piDel); fprintf(fp, "\n");
    
    fclose(fp);
  }
  
  /* parameters */
  
  if(cfg->parms){
    sprintf(file, "parms-chain%d.txt", cfg->chainNum);
    fp = fopen(file, "w");
    
    if(fp == NULL){
      fprintf(stderr, "ERROR: unable to create file, %s\n", file);
      exit(EXIT_FAILURE);
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
      
    for(n = 0; n < cfg->N; ++n)  
      for(g = 0; g < cfg->G; ++g) 
        fprintf(fp, "epsilon_lib%d_gene%d ", n, g);

    fprintf(fp, "\n");   
    
    /* print initial values */

    for(n = 0; n < cfg->N; ++n){
      fprintf(fp, NUM_TF, a->c[n]);
      fprintf(fp, " ");
    }
    
    for(g = 0; g < cfg->G; ++g){
      fprintf(fp, NUM_TF, a->phi[g]);
      fprintf(fp, " ");
    }
    
    for(g = 0; g < cfg->G; ++g){
      fprintf(fp, NUM_TF, a->alp[g]);
      fprintf(fp, " ");
    }
    
    for(g = 0; g < cfg->G; ++g){
      fprintf(fp, NUM_TF, a->del[g]);
      fprintf(fp, " ");
    }
    
    for(g = 0; g < cfg->G; ++g){
      fprintf(fp, NUM_TF, a->eta[g]);
      fprintf(fp, " ");
    }
    
    for(n = 0; n < cfg->N; ++n)
      for(g = 0; g < cfg->G; ++g){  
        fprintf(fp, NUM_TF, a->eps[iG(n, g)]);
        fprintf(fp, " ");
      }
     
    fprintf(fp, "\n");  
    fclose(fp);
  }
  
  /* acceptance rates of Metropolis steps */
  
  if(cfg->rates){
  
    sprintf(file, "rates-chain%d.txt", cfg->chainNum);
    fp = fopen(file, "w"); 
    
    if(fp == NULL){
      fprintf(stderr, "ERROR: unable to create file, %s\n", file);
      exit(EXIT_FAILURE);
    }
       
    fprintf(fp, "d c phi alpha delta mean-epsilon\n");
    fclose(fp);
  }
  
  /* time spent sampling each parameter */
  
  if(cfg->time){
  
    sprintf(file, "time-chain%d.txt", cfg->chainNum);
    fp = fopen(file, "w"); 
  
    if(fp == NULL){
      fprintf(stderr, "ERROR: unable to create file, %s\n", file);
      exit(EXIT_FAILURE);
    }
    
    fprintf(fp, "c tau pi-alpha pi-delta d theta-phi theta-alpha theta-delta sigma-c ");
    fprintf(fp, "sigma-phi sigma-alpha sigma-delta eta epsilon phi alpha delta\n");
    
    fclose(fp);  
  }
}