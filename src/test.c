#include <Chain.h>
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
     sampleD(a, cfg);
     samplePhiAlpDel(a, cfg);
     sampleSigC(a);
     sampleEta(a);
     sampleTau(a, cfg);
     sampleThePhi(a, cfg);
     sampleTheAlp(a, cfg);
     sampleTheDel(a, cfg);
     sampleSigPhi(a, cfg);
     sampleSigAlp(a, cfg);
     sampleSigDel(a, cfg);
     samplePiAlp(a, cfg);
     samplePiDel(a, cfg);
   }
  
  printChain(a);
  
  freeConfig(cfg);
  freeChainHost(a, cfg);
  return 0;
}