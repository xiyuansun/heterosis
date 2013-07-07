#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <numericTypes.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv){
  Config *cfg = config(argc, argv);
  
  Chain *a = newChainHost(cfg);
  if(a == NULL){
    freeConfig(cfg);
    exit(EXIT_FAILURE);
  }
  
  printConfig(cfg);
  
  printf("\n-----\n");
  printChain(a);
  
  freeConfig(cfg);
  freeChainHost(a, cfg);
  return 0;
}