#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

void oneChain(int argc, char **argv){
  Config *cfg = config(argc, argv);  
  Chain *a = newChain(cfg);
  
  printf("hi1\n");
  
  if(a == NULL){
    freeConfig(cfg);
    exit(EXIT_FAILURE);
  }
    printf("hi2\n");
  
  runChain(a, cfg);

  printf("hi3\n");

  summarizeChain(a, cfg);  

  printf("hi4\n");

  freeChain(a, cfg);
  
    printf("hi5\n");
  
  freeConfig(cfg);
  
    printf("hi6\n");
}