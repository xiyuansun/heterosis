#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

void oneChain(int argc, char **argv){

  Config *cfg = config(argc, argv); 
  Chain *host_a, *dev_a;
  
  if(a == NULL){
    freeConfig(cfg);
    exit(EXIT_FAILURE);
  }
  
  newChain(&host_a, &dev_a, cfg); 
  
  runChain(dev_a, cfg);
  summarizeChain(host_a, dev_a, cfg);
  
  printChain(host_a, dev_a, cfg);
  
  freeChain(host_a, cfg, 0); 
  cudaFree(dev_a);
  freeConfig(cfg);
}