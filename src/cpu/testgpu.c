#include <Config.h>
#include <gpuFunctions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv){
  Config *cfg = config(argc, argv);
  printConfig(cfg);
  return 0;
}
