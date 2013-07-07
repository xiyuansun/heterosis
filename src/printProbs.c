#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <Summary.h>
#include <stdio.h>
#include <stlib.h>

void printProbs(Chain *a, Config *cfg){
  int m, n, g, i, nlibs, ngenes, niter = cfg->M - cfg->burnin;
  num_t tmp, phi, alp, del, alp2, delp2;
  num_t prob_de, prob_hph, prob_lph, prob_mph;
  num_t accD, accC, accPhi, accAlp, accDel, acceps;
  FILE *fp;
  
  if(cfg->probsFlag){
    fp = fopen(cfg->probsFile, "w");
    
    fprintf("de ");
    if(cfg->heterosis)
      fprintf("hph lph mph");
    fprintf("\n");
    
    for(g = 0; g < cfg->G; ++g){
      prob_de = 0;
      
      for(m = cfg->burnin + 1; m <= cfg->M; ++m){
        alp = a->alp[m][g];
        prob_de += ((alp * alp) > 1e-6);
      }
      
      prob_de /= niter;
      fprintf(fp, NUM_TF, prob_de);
      fprintf(fp, " ");
      
      if(cfg->heterosis){ 
      
        prob_hph = 0;
        prob_lph = 0;
        prob_mph = 0;
      
        for(m = cfg->burnin + 1; m <= cfg->M; ++m){
        
          phi = a->phi[m][g];
          alp = a->alp[m][g];
          del = a->del[m][g];
          
          alp2 = alp * alp;
          del2 = del * del;
          
          prob_hph += (del2 >  alp2);
          prob_lph += (del2 < -alp2);
          prob_mph += (del2 > 1e-6);
        }
      
        prob_hph /= niter;
        prob_lph /= niter;
        prob_mph /= niter;
        
        fprintf(NUM_TF, prob_hph); fprintf(" ");
        fprintf(NUM_TF, prob_lph); fprintf(" ");
        fprintf(NUM_TF, prob_mph); fprintf(" ");
      }
      fprintf("\n");
    }
    
    fclose(fp);
  }
  
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