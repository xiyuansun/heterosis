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
  Chain *a = allocChain(cfg);
  freeChain(a, cfg);
 
  return EXIT_SUCCESS;
}