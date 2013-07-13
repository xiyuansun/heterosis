#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

__host__ void chains(int argc, char **argv){
  int i;
  Config *cfg = config(argc, argv); 
  
  if(cfg->verbose)  
    printf("Running %d chain(s).\n", cfg->chains);
  
  for(i = 0; i < cfg->chains; ++i){
    if(cfg->verbose) 
      printf("Chain %d\n", i);

    cfg->chainNum = i;
    oneChain(cfg);
  }

  if(cfg->verbose)
    printf("Done running mcmc.\n");
  freeConfig(cfg);
}