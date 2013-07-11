#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

__host__ void chains(int argc, char **argv){
  int i;
  Config *cfg = config(argc, argv); 
  fprintf(cfg->log, "Begun chains");
  
  for(i = 0; i < cfg->M; ++i){
  printf("%d\n", i);
    fprintf(cfg->log, "Chain %d\n", i);
    cfg->chainNum = i;
    oneChain(cfg);
  }

  fprintf(cfg->log, "Done running mcmc.\n");
  free(cfg);
}