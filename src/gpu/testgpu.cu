#include <Config.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv){
  Config *cfg = config(argc, argv);
  Chain *a = allocChain(cfg);
  return 0;
}
