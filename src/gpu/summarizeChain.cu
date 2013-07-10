#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

__host__ void summarizeChain(Chain *host_a, Chain dev_a, Config *cfg){
  Chain *allHost_a = chainDeviceToHost(host_a, dev_a, cfg);

  printProbs(allHost_a, cfg);
  printRates(allHost_a, cfg);
  printHyper(allHost_a, cfg);
  printParms(allHost_a, cfg);
  
  freeChain(allHost_a, cfg, 1);
}