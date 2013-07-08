#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

void oneChain(int argc, char **argv){

  Config *cfg = config(argc, argv);  
  Chain *a = newChain(cfg);
  
  if(a == NULL){
    freeConfig(cfg);
    exit(EXIT_FAILURE);
  }
  
  printf("1\n");
  
  runChain(a, cfg);
  
  printf("2\n");
  
  summarizeChain(a, cfg);  

  freeChain(a, cfg);
  freeConfig(cfg);
}