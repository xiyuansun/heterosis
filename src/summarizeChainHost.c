#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <Summary.h>
#include <stdio.h>
#include <stlib.h>

void summarizeChain(Chain *a, Config *cfg){
  printProbs(a, cfg);
  printRates(a, cfg);
  printHyper(a, cfg);
  printParms(a, cfg);
}