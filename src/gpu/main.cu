#include <Chain.h>
#include <Config.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv){
  
  Config *cfg = config(argc, argv); 
  Chain *host_a, *dev_a;
  
  newChain(&host_a, &dev_a, cfg); 
  printChain(host_a, dev_a, cfg);
  
  freeChain(host_a, cfg, 0); 
  cudaFree(dev_a);
 
  return EXIT_SUCCESS;
}