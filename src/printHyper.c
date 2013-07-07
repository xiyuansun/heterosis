#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <Summary.h>
#include <stdio.h>
#include <stlib.h>

void printHyper(Chain *a, Config *cfg){

  int m, n, g, i, nlibs, ngenes, niter = cfg->M - cfg->burnin;
  num_t tmp, phi, alp, del, alp2, delp2;
  num_t prob_de, prob_hph, prob_lph, prob_mph;
  num_t accD, accC, accPhi, accAlp, accDel, acceps;
  FILE *fp;
  
  if(cfg->hyperFlag){
    fp = fopen(cfg->hyperFile, "w");
    
    fprintf(fp, "sigC d tau thePhi theAlp theDel sigPhi sigAlp sigDel piAlp piDel\n");
    
    for(m = 0; m <= cfg->M; ++m){
      tmp = a->sigC[m];   
      fpritnf(fp, NUM_TF, tmp); fprintf(" "); 

      tmp = a->d[m];      
      fpritnf(fp, NUM_TF, tmp); fprintf(" "); 
      
      tmp = a->tau[m];    
      fpritnf(fp, NUM_TF, tmp); fprintf(" "); 
      
      tmp = a->thePhi[m]; 
      fpritnf(fp, NUM_TF, tmp); fprintf(" "); 
      
      tmp = a->theAlp[m]; 
      fpritnf(fp, NUM_TF, tmp); fprintf(" "); 
      
      tmp = a->theDel[m]; 
      fpritnf(fp, NUM_TF, tmp); fprintf(" "); 
      
      tmp = a->sigPhi[m]; 
      fpritnf(fp, NUM_TF, tmp); fprintf(" "); 
      
      tmp = a->sigAlp[m]; 
      fpritnf(fp, NUM_TF, tmp); fprintf(" "); 
      
      tmp = a->sigDel[m]; 
      fpritnf(fp, NUM_TF, tmp); fprintf(" "); 
      
      tmp = a->piAlp[m];  
      fpritnf(fp, NUM_TF, tmp); fprintf(" "); 
      
      tmp = a->piDel[m];  
      fpritnf(fp, NUM_TF, tmp); fprintf(" "); 
      
      fprintf(fp, "\n");
    }
    
    fclose(fp);
  }
}