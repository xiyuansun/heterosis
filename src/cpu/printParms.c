#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void printParms(Chain *a, Config *cfg){
    
  int m, n, g;
  int N = a->N, G = a->G;
  char file[BUF];
  num_t tmp;
  FILE *fp;
  double time;
  clock_t start = clock();
  
  if(cfg->parmsFlag){   
    fprintf(cfg->log, "  Printing parameters.\n");  
      
    sprintf(file, "../out/parms/chain%d.txt", cfg->chainNum);
    fp = fopen(file, "w");
      
    if(fp == NULL){
      printf("ERROR: unable to create file, %s\n", file);
      return;
    }
    
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
      
    for(g = 0; g < cfg->G; ++g)
      for(n = 0; n < cfg->N; ++n)
        fprintf(fp, "eps_lib%d_gene%d ", n, g);

    fprintf(fp, "\n");
    
    for(m = 0; m <= cfg->M; ++m){
      for(n = 0; n < cfg->N; ++n){
        tmp = a->c[iN(m, n)];
        fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
      }
      
      for(g = 0; g < cfg->G; ++g){
        tmp = a->phi[iG(m, g)];
        fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
      }

      for(g = 0; g < cfg->G; ++g){
        tmp = a->alp[iG(m, g)];
        fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
      }
      
      for(g = 0; g < cfg->G; ++g){
        tmp = a->del[iG(m, g)];
        fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
      }      

      for(g = 0; g < cfg->G; ++g){
        tmp = a->eta[iG(m, g)];
        fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
      }    
      
      for(n = 0; n < cfg->N; ++n)
        for(g = 0; g < cfg->G; ++g){
          tmp = a->eps[iNG(m, n, g)];
          fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
        }
      
      fprintf(fp, "\n");
    } 
    
    fclose(fp);
  }
  
  time = ((double) clock() - start) / (60 * CLOCKS_PER_SEC);
  fprintf(cfg->time, "%0.3f ", time);
}