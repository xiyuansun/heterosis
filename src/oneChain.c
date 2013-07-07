#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

void oneChain(int argc, char **argv){
  Config *cfg = config(argc, argv);  
  Chain *a = newChainHost(cfg);
  
  if(a == NULL){
    freeConfig(cfg);
    exit(EXIT_FAILURE);
  }
  
  runChain(a, cfg);
  summarizeChain(a, cfg);
  printChain(a);

  freeConfig(cfg);
  freeChainHost(a, cfg);
}