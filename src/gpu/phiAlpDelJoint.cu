#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

__device__ num_t lPhiAlpDelJoint(Chain *a, int g, num_t argPhi, num_t argAlp, num_t argDel){ /* device */
 
  int n, G = a->G;
  num_t ret, s = 0, tmp = 0;

  for(n = 0; n < a->N; ++n){
    tmp = mu(a, n, argPhi, argAlp, argDel);
    s += a->y[iG(n, g)] * tmp - exp(a->c[n] + 
         a->eps[iG(n, g)] + tmp);
  }

  /* phi part */
  ret = s - pow(argPhi - a->thePhi, 2) / (2 * pow(a->sigPhi, 2));

  /* alpha part */
  if(argAlp * argAlp > 1e-6){
    tmp = -pow(argAlp - a->theAlp, 2) / (2 * pow(a->sigAlp, 2)) -
                log(1 - a->piAlp);
  } else {
    tmp = log(a->piAlp);
  }

  ret = ret + tmp;

  /* delta part */
  if(argDel * argDel > 1e-6){
    tmp = -pow(argDel - a->theDel, 2) / (2 * pow(a->sigDel, 2)) -
                log(1 - a->piDel);
  } else {
    tmp = log(a->piDel);
  }

  return ret + tmp;
}

__global__ void samplePhiAlpDelJoint_kernel(Chain *a){ /* kernel <<<G, 1>>> */
  int g = IDX;
  num_t oldPhi, newPhi, oldAlp, newAlp, oldDel, newDel;
  num_t dl, lp, lu;

  if(g < G){

    oldPhi = a->phi[g];
    newPhi = rnormalDevice(a, g, oldPhi, a->tunePhi[g]);

    oldAlp = a->alp[g];
    newAlp = alpProp(a, g);

    oldDel = a->del[g];
    newDel = delProp(a, g);

    dl = lPhiAlpDelJoint(a, g, newPhi, newAlp, newDel) 
       - lPhiAlpDelJoint(a, g, oldPhi, oldAlp, oldDel); 
    lp = 0 < dl ? 0 : dl;
    lu = log(runiformDevice(a, g, 0, 1));

    if(lu < lp){ /* accept */
      a->phi[g] = newPhi;
      a->alp[g] = newAlp;
      a->del[g] = newDel;

      a->tunePhi[g] *= 1.1; 

      if(a->m > a->burnin){
        ++a->accPhi[g];
        ++a->accAlp[g];
        ++a->accDel[g];
      }
    } else { /* reject */
      a->tunePhi[g] /= 1.1;
    }
  }
}

__host__ void samplePhiAlpDelJoint(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */

  float myTime;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  if(cfg->verbose)
    printf("phiAlpDelJoint ");

  samplePhiAlpDelJoint_kernel<<<G_GRID, G_BLOCK>>>(dev_a);

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  
  cfg->timePhi = myTime;
  cfg->timeAlp = myTime;
  cfg->timeDel = myTime;  
}