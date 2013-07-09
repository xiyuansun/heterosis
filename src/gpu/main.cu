#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv){
 /* oneChain(argc, argv);*/
 
  Config *cfg = config(argc, argv);  
  Chain *a = allocChain(cfg);
  freeChain(a);
 
  return EXIT_SUCCESS;
}