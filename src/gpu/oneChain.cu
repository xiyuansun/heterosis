#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

void oneChain(Config*){

  Chain *host_a = NULL, *dev_a = NULL;
  
  newChain(&host_a, &dev_a, cfg); 
  
  if(host_a == NULL || dev_a == NULL){
    freeConfig(cfg);
    exit(EXIT_FAILURE);
  }
  
  runChain(host_a, dev_a, cfg);
  summarizeChain(host_a, dev_a, cfg);
  
  freeChain(host_a, cfg, 0); 
  cudaFree(dev_a);
}