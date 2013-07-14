#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

void resetChain(Chain *a, Config *cfg){
  copyHyperParms(a, cfg);
  newChain_kernel1(a);
  newChain_kernel2(a);
}