#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

void samplePhiAlpDel(Chain *a, Config *cfg){ /* host */
  if(cfg->joint && cfg->heterosis){
    samplePhiAlpDelJoint(a);
  } else {
    samplePhi(a);
    sampleAlp(a);
    
    if(cfg->heterosis)
      sampleDel(a);
  }
} 