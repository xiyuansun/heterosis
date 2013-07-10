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
 
  Config *cfg = config(argc, argv); 
  Chain *a = newChain(cfg); 
/*  printChain(a, cfg, 0); */
 /* freeChain(a, cfg, 0); */
 
  return EXIT_SUCCESS;
}