#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

void breakpoint(){
  return;
}

void samplePhiAlpDel(Chain *a, Config *cfg){ /* host */


if(a->m == 6){
  breakpoint();
}

  if(cfg->joint && cfg->heterosis){
    samplePhiAlpDelJoint(a, cfg);
  } else {
    samplePhi(a, cfg); 
    sampleAlp(a, cfg); 
    
    if(cfg->heterosis)
      sampleDel(a, cfg);
  }
} 