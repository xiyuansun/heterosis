#include <Chain.h>
#include <Config.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv){
 /* oneChain(argc, argv);*/
 
  Config *cfg = config(argc, argv);  printf("1\n");
  Chain *a = allocChain(cfg); printf("2\n");
  printChain(a, cfg, 0); printf("3\n");
  freeChain(a, cfg, 0);printf("4\n");
 
  return EXIT_SUCCESS;
}