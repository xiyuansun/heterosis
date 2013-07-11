#include <Config.h>
#include <functions.h>
#include <stdlib.h>

void freeConfig(Config *cfg){
  free(cfg);
}