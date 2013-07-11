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
  float myTime;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  
  fprintf(cfg->log, "  Running chain.\n");
  
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
  
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  
  fprintf(cfg->time, "%0.3f ", myTime/60000.0); /* elapsed time in minutes */
} 