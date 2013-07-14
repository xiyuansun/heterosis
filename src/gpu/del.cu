#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

__device__ num_t lDel(Chain *a, int g, num_t arg){ /* device */ 
  int n, G = a->G;
  num_t s = 0, tmp; 
  
  for(n = 0; n < a->N; ++n){
    if(a->grp[n] != 2){
      tmp = mu(a, n, a->phi[g], a->alp[g], arg);
      s += a->y[iG(n, g)] * tmp - exp(a->c[n] + 
          a->eps[iG(n, g)] + tmp);
    }
  }
 
  if(arg * arg > 1e-6){
    tmp = -pow((float) (arg - a->theDel), 2) / (2 * pow((float) a->sigDel, 2)) -
                log(1 - a->piDel);
  } else {
    tmp = log(a->piDel);
  }

  return s + tmp;
}

__global__ void sampleDel_kernel(Chain *a){ /* kernel <<<G, 1>>> */
  int g = IDX;
  num_t old, nw, dl, lp, lu;

  if(g < a->G){ 

    old = a->del[g];
    nw = delProp(a, g);
    
    dl = lDel(a, g, nw) - lDel(a, g, old);
    lp = 0 < dl? 0 : dl;
    lu = log(runiformDevice(a, g, 0, 1));
    
    if(lu < lp){ /* accept */
      a->del[g] = nw;
      
      if(a->m > a->burnin)
        ++a->accDel[g];
    } else { /* reject */
      a->del[g] = old;
    }
  }
}

__host__ void sampleDel(Chain *host_a, Chain *dev_a, Config* cfg){ /* host */

  num_t myTime;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  if(cfg->verbose)
    printf("del ");

  sampleDel_kernel<<<G_GRID, G_BLOCK>>>(dev_a);

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  cfg->timeDel = myTime / MILLISECS;
}