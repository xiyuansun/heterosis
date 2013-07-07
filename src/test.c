#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv){
  int m;

  Config *cfg = config(argc, argv);
  
  Chain *a = newChainHost(cfg);
  if(a == NULL){
    freeConfig(cfg);
    exit(EXIT_FAILURE);
  }
  
 
  printChain(a);
  
   printf("\n-----\n");
   
   for(m = 0; m < a->M; ++m){
     sampleC(a);
     sampleEps(a);
     sampleD(a);
     samplePhiAlpDel(a);
     sampleSigC(a);
     sampleEta(a);
     sampleTau(a);
   }
  
  printChain(a);
  
  freeConfig(cfg);
  freeChainHost(a, cfg);
  return 0;
}