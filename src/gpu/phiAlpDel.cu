#include <Chain.h>
#include <constants.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

void samplePhiAlpDel(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */
  if(cfg->joint && cfg->heterosis){
    samplePhiAlpDelJoint(host_a, dev_a, cfg);
  } else {
    samplePhi(host_a, dev_a, cfg);
    sampleAlp(host_a, dev_a, cfg);
    
    if(cfg->heterosis)
      sampleDel(host_a, dev_a, cfg);
  }
} 