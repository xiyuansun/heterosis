#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

__device__ num_t lEps(Chain *a, int n, int g, num_t arg){ /* device */
  int G = a->G;
  return a->y[iG(n, g)] * arg - exp(a->c[n] + arg + mu(a, n, a->phi[g], 
                                     a->alp[g], a->del[g])) 
                          - (arg * arg) / (2 * pow((float) a->eta[g], 2));
}

__global__ void sampleEps_kernel(Chain *a){ /* kernel <<<N, G>>> */
  int n, g = IDX, G = a->G;
  num_t old, nw, dl, lp, lu;

  if(g < G){
    for(n = 0; n < a->N; ++n){ 
      old = a->eps[iG(n, g)];
      nw = rnormalDevice(a, g, old, a->tuneEps[iG(n, g)]);

      dl = lEps(a, n, g, nw) - lEps(a, n, g, old);
      lp = 0 < dl ? 0 : dl;
      lu = log(runiformDevice(a, g, 0, 1));
      
      if(lu < lp){ /* accept */
        a->eps[iG(n, g)] = nw;
        a->tuneEps[iG(n, g)] *= 1.1;
        
      if(a->m > a->burnin)
          ++a->accEps[iG(n, g)]; 
      } else { /* reject */
        a->eps[iG(n, g)] = old;
        a->tuneEps[iG(n, g)] /= 1.1;
      }
      
      if(a->m > a->burnin)
        a->sumEps[iG(n, g)] += a->eps[iG(n, g)];
    }
  }
}

__host__ void sampleEps(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */

  float myTime;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  if(cfg->verbose)
    printf("eps ");

  sampleEps_kernel<<<G_GRID, G_BLOCK>>>(dev_a);

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  cfg->timeEps = myTime / MILLISECS;
}