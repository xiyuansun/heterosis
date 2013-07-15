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

  newChain(&host_a, &dev_a, cfg);  
  
  if(host_a == NULL){
    free(cfg);
    exit(EXIT_FAILURE);
  }
  
  if(cfg->verbose)
    printf("Running %d chain(s).\n", cfg->chains);
  
  printConfig(cfg);
  
  for(i = 0; i < cfg->chains; ++i){
  
  printChain(host_a, dev_a, cfg);
  
    if(cfg->verbose)
      printf("  Chain %d\n", i);

    cfg->chainNum = i;
    
    runChain(host_a, dev_a, cfg);
    summarizeChain(host_a, dev_a, cfg);

    resetChain(host_a, dev_a, cfg);
  }
  
  freeChain(host_a, cfg, 0);
  cudaFree(dev_a);
  
  if(cfg->gelman)
    gelmanFactors(cfg);
  
  if(cfg->verbose)
    printf("Done.\n");

  free(cfg);
}