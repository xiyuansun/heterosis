#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

void mcmc(int argc, char **argv){
  int i;
  Config *cfg = config(argc, argv); 
  Chain *a = newChain(cfg);
  
  if(a == NULL){
    free(cfg);
    exit(EXIT_FAILURE);
  }
  
  if(cfg->verbose)
    printf("Running %d chain(s).\n", cfg->chains);
  
  for(i = 0; i < cfg->chains; ++i){

    if(cfg->verbose)
      printf("Chain %d\n", i);

    cfg->chainNum = i;
    
    runChain(a, cfg);
    summarizeChain(a, cfg);
    
    resetChain(a, cfg);
  }
  
  freeChain(a, cfg);
  
  if(cfg->verbose)
    printf("Done.\n");
    
  free(cfg);
}