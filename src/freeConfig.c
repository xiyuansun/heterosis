#include <Config.h>
#include <functions.h>
#include <stdlib.h>

void freeConfig(Config *cfg){
  free(cfg->dataFile);
  free(cfg->groupFile);
  free(cfg->probsFile);
  free(cfg->hyperFile);
  free(cfg->ratesFile);
  free(cfg->parmsFile);
  free(cfg);
}