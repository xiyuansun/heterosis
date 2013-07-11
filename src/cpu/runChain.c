#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void runChain(Chain *a, Config *cfg){
  int m;
  fprintf(cfg->log, "  Running chain.\n");
  
  for(m = 0; m < a->M; ++m){
    fprintf(cfg->log, "    iter %d | ", m);
    
    sampleC(a, cfg);
    sampleTau(a, cfg);
    samplePiAlp(a, cfg);
    samplePiDel(a, cfg);
    sampleD(a, cfg);
    sampleThePhi(a, cfg);
    sampleTheAlp(a, cfg);
    sampleTheDel(a, cfg);
    sampleSigC(a, cfg);
    sampleSigPhi(a, cfg);
    sampleSigAlp(a, cfg);
    sampleSigDel(a, cfg);
    sampleEta(a, cfg);
    sampleEps(a, cfg);
    samplePhiAlpDel(a, cfg);
    
    fprintf(cfg->log, "\n");
  }
}