#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

__host__ void oneIteration(Chain *host_a, Chain *dev_a, Config *cfg){
  ++cfg->chainNum;

  if(cfg->verbose)
    printf("  Chain %d of %d.\n", cfg->chainNum, cfg->chains);
      
  runChain(host_a, dev_a, cfg);
  summarizeChain(host_a, dev_a, cfg);
  resetChain(host_a, dev_a, cfg);
}

void mcmc(int argc, char **argv){
  int i;
  Config *cfg = config(argc, argv); 
  Chain *host_a = NULL, *dev_a = NULL;

printConfig(cfg);

  newChain(&host_a, &dev_a, cfg);  
  
  if(host_a == NULL){
    free(cfg);
    exit(EXIT_FAILURE);
  }
  
  if(cfg->verbose)
    printf("Running %d chain(s).\n", cfg->chains);
  
  for(i = 0; i < cfg->chains; ++i)
    oneIteration(host_a, dev_a, cfg);
  
  freeChain(host_a, cfg, 0);
  cudaFree(dev_a);
  
  if(cfg->verbose)
    printf("Done.\n");

  free(cfg);
} 