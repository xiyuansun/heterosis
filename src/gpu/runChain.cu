#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void runChain(Chain *host_a, Chain *dev_a, Config *cfg){
  int m;
  fprintf(cfg->log, "  Running chain.\n"); 
  
  for(m = 0; m < cfg->M; ++m){

    fprintf(cfg->log, "    iter %d | ", m);

    sampleC(host_a, dev_a, cfg); printf("c\n");
    sampleTau(host_a, dev_a, cfg);printf("tau\n");
    samplePiAlp(host_a, dev_a, cfg);printf("pialp\n");
    samplePiDel(host_a, dev_a, cfg); printf("pidel\n");
    sampleD(host_a, dev_a, cfg); printf("d\n");
    sampleThePhi(host_a, dev_a, cfg);printf("thephi\n");
    sampleTheAlp(host_a, dev_a, cfg);printf("thealp\n");
    sampleTheDel(host_a, dev_a, cfg);printf("thedel\n");
    sampleSigC(host_a, dev_a, cfg);printf("sigc\n");
    sampleSigPhi(host_a, dev_a, cfg);printf("sigphi\n");
    sampleSigAlp(host_a, dev_a, cfg);printf("sigalp\n");
    sampleSigDel(host_a, dev_a, cfg);printf("sigdel\n");
    sampleEta(host_a, dev_a, cfg);  printf("eta\n");
    sampleEps(host_a, dev_a, cfg);printf("eps\n");
    samplePhiAlpDel(host_a, dev_a, cfg);printf("phialpdel\n");

    fprintf(cfg->time, "\n");
    fprintf(cfg->log, "\n");
  }
} 