#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
 
void runChain(Chain *host_a, Chain *dev_a, Config *cfg){
  int m;
  
  if(cfg->verbose)
    printf("  Running chain.\n");
  
  printHeaders(host_a, dev_a, cfg);
  
  for(m = 0; m < cfg->M; ++m){
    
    if(cfg->verbose)
      printf("    iter %d | ", m);
    
    sampleC(host_a, dev_a,  cfg);
    sampleTau(host_a, dev_a,  cfg);
    samplePiAlp(host_a, dev_a, cfg);
    samplePiDel(host_a, dev_a, cfg);
    sampleD(host_a, dev_a,  cfg); 
    sampleThePhi(host_a, dev_a,  cfg);
    sampleTheAlp(host_a, dev_a,  cfg);
    sampleTheDel(host_a, dev_a,  cfg);
    sampleSigC(host_a, dev_a,  cfg);
    sampleSigPhi(host_a, dev_a,  cfg);
    sampleSigAlp(host_a, dev_a, cfg);
    sampleSigDel(host_a, dev_a,  cfg);
    sampleEta(host_a, dev_a,  cfg);
    sampleEps(host_a, dev_a,  cfg);
    samplePhiAlpDel(host_a, dev_a, cfg); 
    
    interimResults(host_a, dev_a, cfg);
    
    if(cfg->verbose)
      printf("\n");
  }
}