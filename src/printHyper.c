#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <Summary.h>
#include <stdio.h>
#include <stdlib.h>

void printHyper(Chain *a, Config *cfg){

  int m;
  num_t tmp;
  FILE *fp;
  
  if(cfg->hyperFlag){
    fp = fopen(cfg->hyperFile, "w");
    
    fprintf(fp, "sigC d tau thePhi theAlp theDel sigPhi sigAlp sigDel piAlp piDel\n");
    
    for(m = 0; m <= cfg->M; ++m){
      tmp = a->sigC[m];   
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " "); 

      tmp = a->d[m];      
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " "); 
      
      tmp = a->tau[m];    
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " "); 
      
      tmp = a->thePhi[m]; 
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " "); 
      
      tmp = a->theAlp[m]; 
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " "); 
      
      tmp = a->theDel[m]; 
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " "); 
      
      tmp = a->sigPhi[m]; 
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " "); 
      
      tmp = a->sigAlp[m]; 
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " "); 
      
      tmp = a->sigDel[m]; 
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " "); 
      
      tmp = a->piAlp[m];  
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " "); 
      
      tmp = a->piDel[m];  
      fprintf(fp, NUM_TF, tmp); fprintf(fp, " "); 
      
      fprintf(fp, "\n");
    }
    
    fclose(fp);
  }
}