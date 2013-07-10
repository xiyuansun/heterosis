#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void runChain(Chain *host_a, Chain *dev_a, Config *cfg){

  int m;
  
  for(m = 0; m < cfg->M; ++m){
    sampleC(host_a, dev_a, cfg);
  /*  sampleTau(a, cfg);
    samplePiAlp(a, cfg);
    samplePiDel(a, cfg);
    sampleD(a, cfg);
    sampleThePhi(a, cfg);
    sampleTheAlp(a, cfg);
    sampleTheDel(a, cfg);
    sampleSigC(a);
    sampleSigPhi(a, cfg);
    sampleSigAlp(a, cfg);
    sampleSigDel(a, cfg);
    sampleEta(a);
    sampleEps(a);
    samplePhiAlpDel(a, cfg);*/
  }
}