#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void runChain(Chain *host_a, Chain *dev_a, Config *cfg){
  int m;
  double time;
  time_t p1, p2;
  
  fprintf(cfg->log, "  Running chain.\n");
  time(p1);
  
  for(m = 0; m < cfg->M; ++m){
    fprintf(cfg->log, "    iter %d | ", m);

    sampleC(host_a, dev_a, cfg);
    sampleTau(host_a, dev_a, cfg);
    samplePiAlp(host_a, dev_a, cfg);
    samplePiDel(host_a, dev_a, cfg);
    sampleD(host_a, dev_a, cfg);
    sampleThePhi(host_a, dev_a, cfg);
    sampleTheAlp(host_a, dev_a, cfg);
    sampleTheDel(host_a, dev_a, cfg);
    sampleSigC(host_a, dev_a, cfg);
    sampleSigPhi(host_a, dev_a, cfg);
    sampleSigAlp(host_a, dev_a, cfg);
    sampleSigDel(host_a, dev_a, cfg);
    sampleEta(host_a, dev_a, cfg);
    sampleEps(host_a, dev_a, cfg);
    samplePhiAlpDel(host_a, dev_a, cfg);

    fprintf(cfg->log, "\n");
  }
  
  time(p2);
  fprintf(cfg->time, "%0.3f ", difftime(p2, p1)/60.0);
} 