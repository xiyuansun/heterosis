#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

__device__ num_t lAlp(Chain *a, int g, num_t arg){ /* device */
  
  int n, G = a->G;
  num_t s = 0, tmp;
   
  for(n = 0; n < a->N; ++n){
    if(a->grp[n] != 2){
      tmp = mu(a, n, a->phi[g], arg, a->del[g]);
      s += a->y[iG(n, g)] * tmp - exp(a->c[n] + 
          a->eps[iG(n, g)] + tmp);
    }
  }
 
  if(arg * arg > 1e-6){
    tmp = -pow((float) (arg - a->theAlp), 2) / (2 * pow((float) a->sigAlp, 2)) -
                log(1 - a->piAlp);
  } else {
    tmp = log(a->piAlp);
  }

  return s + tmp;
}

__global__ void sampleAlp_kernel(Chain *a){ /* kernel <<<G, 1>>> */

  int g = IDX;
  num_t old, nw, dl, lp, lu;

  if(g < a->G){ 

    old = a->alp[g];
    nw = alpProp(a, g);
    
    dl = lAlp(a, g, nw) - lAlp(a, g, old);
    lp = 0 < dl ? 0 : dl;
    lu = log(runiformDevice(a, g, 0, 1));
    
    if(lu < lp){ /* accept */
      a->alp[g] = nw;
      
      if(a->m > a->burnin)
        ++a->accAlp[g];
    } 
  }
}

__host__ void sampleAlp(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */

  num_t myTime;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  if(cfg->verbose)
    printf("alpha ");

  sampleAlp_kernel<<<G_GRID, G_BLOCK>>>(dev_a);
  
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  
  cfg->timeAlp = myTime;
}