#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

void mcmc(int argc, char **argv){
  int i;
  Config *cfg = config(argc, argv); 
  Chain *host_a = NULL, *dev_a = NULL;
  
  if(cfg->verbose)
    printf("Running %d chain(s).\n", cfg->chains);

  newChain(&host_a, &dev_a, cfg);  
  
  if(host_a == NULL){
    free(cfg);
    exit(EXIT_FAILURE);
  }
  
  for(i = 0; i < cfg->chains; ++i){
  
    if(cfg->verbose)
      printf("  Chain %d\n", i);

    cfg->chainNum = i;
    
    printChain(host_a, dev_a, cfg);
    
    runChain(host_a, dev_a, cfg);
    summarizeChain(host_a, dev_a, cfg);

    resetChain(host_a, dev_a, cfg);
  }
  
  freeChain(host_a, cfg, 0);
  cudaFree(dev_a);
  
  if(cfg->verbose)
    printf("Done running mcmc.\n");

  free(cfg);
}