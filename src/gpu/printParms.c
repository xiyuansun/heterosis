#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

void printParms_oneFile(Chain *a, Config *cfg, int some){

  int m, n, g, nlibs, ngenes;
  int N = a->N, G = a->G;
  num_t tmp;
  FILE *fp;
  
  if(cfg->someParmsFlag || cfg->allParmsFlag){   
    if(cfg->someParmsFlag && some){
      fp = fopen(cfg->someParmsFile, "w");
      
      if(fp == NULL){
        printf("ERROR: unable to create file, %s\n", cfg->someParmsFile);
        return;
      }
      
      nlibs = 5 < cfg->N ? 5 : cfg->N;
      ngenes = 5 < cfg->G ? 5 : cfg->G;
    } else if(cfg->allParmsFlag && !some){
      fp = fopen(cfg->allParmsFile, "w");
      
      if(fp == NULL){
        printf("ERROR: unable to create file, %s\n", cfg->allParmsFile);
        return;
      }
      
      nlibs = cfg->N;
      ngenes = cfg->G;
    } else {
      return;
    }
    
    for(n = 0; n < nlibs; ++n)
      fprintf(fp, "c%d ", n);
    
    for(g = 0; g < ngenes; ++g)
      fprintf(fp, "phi%d ", g);
    
    for(g = 0; g < ngenes; ++g)
      fprintf(fp, "alpha%d ", g);
    
    for(g = 0; g < ngenes; ++g)
      fprintf(fp, "delta%d ", g);
      
    for(g = 0; g < ngenes; ++g)
      fprintf(fp, "eta%d ", g);
      
    for(g = 0; g < ngenes; ++g)
      for(n = 0; n < nlibs; ++n)
        fprintf(fp, "eps_lib%d_gene%d ", n, g);

    fprintf(fp, "\n");
    
    for(m = 0; m <= cfg->M; ++m){
      for(n = 0; n < nlibs; ++n){
        tmp = a->c[iN(m, n)];
        fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
      }
      
      for(g = 0; g < ngenes; ++g){
        tmp = a->phi[iG(m, g)];
        fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
      }

      for(g = 0; g < ngenes; ++g){
        tmp = a->alp[m][g];
        fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
      }
      
      for(g = 0; g < ngenes; ++g){
        tmp = a->del[m][g];
        fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
      }      

      for(g = 0; g < ngenes; ++g){
        tmp = a->eta[iG(m, g)];
        fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
      }    
      
      for(n = 0; n < nlibs; ++n)
        for(g = 0; g < ngenes; ++g){
          tmp = a->eps[iNG(m, n, g)];
          fprintf(fp, NUM_TF, tmp); fprintf(fp, " ");
        }
      
      fprintf(fp, "\n");
    } 
    
    fclose(fp);
  }
}

void printParms(Chain *a, Config *cfg){
  printParms_oneFile(a, cfg, 0);
  printParms_oneFile(a, cfg, 1);
}