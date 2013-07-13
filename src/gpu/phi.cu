#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

__device__ num_t lPhi(Chain *a, int g, num_t arg){ /* device */
  int n, G = a->G;
  num_t ret, s = 0, tmp = 0; 

  for(n = 0; n < a->N; ++n){
    tmp = mu(a, n, arg, a->alp[g], a->del[g]);
    s += a->y[iG(n, g)] * tmp - exp(a->c[n] + 
        a->eps[iG(n, g)] + tmp);
  }
 
  ret = s - pow(arg - a->thePhi, 2) / (2 * pow(a->sigPhi, 2));
  return ret;
}

__global__ void samplePhi_kernel(Chain *a){ /* kernel <<<G, 1>>> */
  int g = IDX;
  num_t old, nw, dl, lp, lu;
  
  if(g < a->G){ 

    old = a->phi[g];
    nw = rnormalDevice(a, g, old, a->tunePhi[g]);

    dl = lPhi(a, g, nw) - lPhi(a, g, old);
    lp = 0 < dl ? 0 : dl;
    lu = log(runiformDevice(a, g, 0, 1)); 
    
    if(lu < lp){ /* accept */
      a->phi[g] = nw;
      a->tunePhi[g] *= 1.1; 
      
      if(a->m > a->burnin)
        ++a->accPhi[g];
    } else { /* reject */
      a->phi[g] = old;
      a->tunePhi[g] /= 1.1; 
    }
  }
}

__host__ void samplePhi(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */

  num_t myTime;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  if(cfg->verbose)
    printf("phi ");

  samplePhi_kernel<<<G_GRID, G_BLOCK>>>(dev_a);

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  cfg->timePhi = myTime;
}
