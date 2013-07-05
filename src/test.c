#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv){
  
  Config *cfg = config(argc, argv);
  printConfig(cfg);
  freeConfig(cfg);

  return 0;
}