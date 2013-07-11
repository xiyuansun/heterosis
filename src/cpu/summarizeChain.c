#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

void summarizeChain(Chain *a, Config *cfg){

  printProbs(a, cfg);
  printRates(a, cfg);
  printHyper(a, cfg);
  printParms(a, cfg);
  
  fprintf(cfg->time, "\n");
}