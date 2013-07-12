#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

__host__ void chains(int argc, char **argv){
  int i;
  Config *cfg = config(argc, argv); 
  fprintf(cfg->log, "Running %d chains.\n", cfg->chains);
  
  for(i = 0; i < cfg->chains; ++i){
    fprintf(cfg->log, "Chain %d\n", i);
    cfg->chainNum = i;
    oneChain(cfg);
  }

  fprintf(cfg->log, "Done running mcmc.\n");
  free(cfg);
}