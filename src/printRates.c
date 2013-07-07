#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <Summary.h>
#include <stdio.h>
#include <stlib.h>

void printRates(Chain *a, Config *cfg){

  int m, n, g, i, nlibs, ngenes, niter = cfg->M - cfg->burnin;
  num_t tmp, phi, alp, del, alp2, delp2;
  num_t prob_de, prob_hph, prob_lph, prob_mph;
  num_t accD, accC, accPhi, accAlp, accDel, acceps;
  FILE *fp;

  if(cfg->ratesFlag){
    fp = fopen(cfg->ratesFile, "w");    
    fprintf(fp, "d c phi alp del meanEps\n");
    
    accD   = a->accD;
    accC   = a->accC[0];
    accPhi = a->accPhi[0];
    accAlp = a->accAlp[0];
    accDel = a->accDel[0];
    accEps = a->accEps[0];
    
    fprintf(fp, NUM_TF, accD);   fprintf(fp, " ");
    fprintf(fp, NUM_TF, accC);   fprintf(fp, " ");
    fprintf(fp, NUM_TF, accPhi); fprintf(fp, " ");
    fprintf(fp, NUM_TF, accAlp); fprintf(fp, " ");
    fprintf(fp, NUM_TF, accDel); fprintf(fp, " ");
    fprintf(fp, NUM_TF, accEps); fprintf(fp, " ");

    for(i = 0; i < cfg->N; ++i){
      accC   = a->accC[i];
      accPhi = a->accPhi[i];
      accAlp = a->accAlp[i];
      accDel = a->accDel[i];
      accEps = a->accEps[i];
      
      fprintf(fp, fprintf(fp, ". ");
      fprintf(fp, NUM_TF, accC);   fprintf(fp, " ");
      fprintf(fp, NUM_TF, accPhi); fprintf(fp, " ");
      fprintf(fp, NUM_TF, accAlp); fprintf(fp, " ");
      fprintf(fp, NUM_TF, accDel); fprintf(fp, " ");
      fprintf(fp, NUM_TF, accEps); fprintf(fp, " ");, 
    }
    
    for(i = cfg->N; i < cfg->G; ++i){
      accPhi = a->accPhi[i];
      accAlp = a->accAlp[i];
      accDel = a->accDel[i];
      accEps = a->accEps[i];
      
      fprintf(fp, fprintf(fp, ". . ");
      fprintf(fp, NUM_TF, accPhi); fprintf(fp, " ");
      fprintf(fp, NUM_TF, accAlp); fprintf(fp, " ");
      fprintf(fp, NUM_TF, accDel); fprintf(fp, " ");
      fprintf(fp, NUM_TF, accEps); fprintf(fp, " ");
    }
    
    fclose(fp);
  }
}