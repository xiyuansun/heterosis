#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void samplePhiAlpDel(Chain *a, Config *cfg){ /* host */
  if(cfg->joint){
    samplePhiAlpDelJoint(a);
  } else {
    samplePhi(a);
    sampleAlp(a);
    sampleDel(a);
  }
} 