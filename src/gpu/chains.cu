#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

__host__ void chains(int argc, char **argv){
  int i;
  Config *cfg = config(argc, argv); 
  printf("okay");
  for(i = 0; i < cfg->M; ++i){
    cfg->chainNum = i;
    oneChain(cfg);
  }
  
  free(cfg);
}