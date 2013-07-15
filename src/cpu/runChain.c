#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void runChain(Chain *a, Config *cfg){
  int m;
  
  if(cfg->verbose)
    printf("  Running chain.\n");
  
  printHeaders(a, cfg);
  
  for(m = 0; m < a->M; ++m){
    
    if(cfg->verbose)
      printf("    iter %d | ", m + 1);
    
    sampleC(a, cfg);
    sampleTau(a, cfg);
    samplePiAlp(a, cfg);
    samplePiDel(a, cfg);
    sampleD(a, cfg); 
    sampleThePhi(a, cfg);
    sampleTheAlp(a, cfg);
    sampleTheDel(a, cfg);
    sampleSigC(a, cfg);
    sampleSigPhi(a, cfg);
    sampleSigAlp(a, cfg);
    sampleSigDel(a, cfg);
    sampleEta(a, cfg);
    sampleEps(a, cfg);
    samplePhiAlpDel(a, cfg); 
    
    intermResults(a, cfg);
    
    if(cfg->verbose)
      printf("\n");
  }
}