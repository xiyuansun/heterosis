#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

void chains(int argc, char **argv){
  Config *cfg = config(argc, argv); 
  
  for(cfg->chainNum = 0; cfg->chainNum < cfg->M; ++cfg->chainNum)
    oneChain(cfg);
  
  free(cfg);
}