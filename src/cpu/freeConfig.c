#include <Chain.h>
#include <Config.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

void freeConfig(Config *cfg){
  fclose(cfg->log);
  free(cfg);
}