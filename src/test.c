#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <numericTypes.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv){
  Config *cfg = config(argc, argv);
  Chain *a = newChainHost(cfg);
  printChain(a);
  
  freeConfig(cfg);
  freeChainHost(a, cfg);
  return 0;
}