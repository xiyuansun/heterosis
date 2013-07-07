#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <Summary.h>
#include <stdio.h>
#include <stlib.h>

void printParms(Chain *a, Config *cfg){

  int m, n, g, i, nlibs, ngenes, niter = cfg->M - cfg->burnin;
  num_t tmp, phi, alp, del, alp2, delp2;
  num_t prob_de, prob_hph, prob_lph, prob_mph;
  num_t accD, accC, accPhi, accAlp, accDel, acceps;
  FILE *fp;
  
  if(cfg->someParmsFlag || allParmsFlag){
    fp = fopen(cfg->hyperFile, "w");
    
    if(cfg->someParmsFlag){
      nlibs = 5;
      ngenes = 5;
    } else {
      nlibs = cfg->N;
      ngenes = cfg->G;
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
        
    fprintf("\n");
    
    for(m = 0; m <= cfg->M; ++m){
      for(n = 0; n < nlibs; ++n){
        tmp = a->c[m][n];
        fprintf(fp, NUM_TF, tmp); fprintf(fp " ");
      }
      
      for(g = 0; g < ngenes; ++g){
        tmp = a->phi[m][g];
        fprintf(fp, NUM_TF, g); fprintf(fp " ");
      }

      for(g = 0; g < ngenes; ++g){
        tmp = a->alp[m][g];
        fprintf(fp, NUM_TF, g); fprintf(fp " ");
      }
      
      for(g = 0; g < ngenes; ++g){
        tmp = a->del[m][g];
        fprintf(fp, NUM_TF, g); fprintf(fp " ");
      }      

      for(g = 0; g < ngenes; ++g){
        tmp = a->eta[m][g];
        fprintf(fp, NUM_TF, g); fprintf(fp " ");
      }      
      
      for(n = 0; n < nlibs; ++n)
        for(g = 0; g < ngenes; ++g){
          tmp = a->eps[m][n][g];
          fprintf(fp, NUM_TF, g); fprintf(fp " ");
        }
      
      fprintf("\n");
    } 
    
    fclose(fp);
  }
}