#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

__host__ void summarizeChain(Chain *host_a, Chain *dev_a, Config *cfg){
  /*printProbs(host_a, cfg);
  printRates(host_a, dev_a, cfg);
  printHyper(host_a, cfg);
  printParms(host_a, dev_a, cfg);*/
  
  fprintf(cfg->timeConfig, "\n");
}