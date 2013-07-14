#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>

__host__ void resetChain(Chain *host_a, Chain *dev_a, Config *cfg){
  newChain_kernel1<<<G_GRID, G_BLOCK>>>(dev_a);
  newChain_kernel2<<<1, 1>>>(dev_a);
}