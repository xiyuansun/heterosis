#include <Chain.h>
#include <constants.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

__device__ num_t lPhiAlpDelJoint(Chain *a, int g, num_t argPhi, num_t argAlp, num_t argDel){ /* device */
 
  int n, N = a->N, G = a->G;
  num_t ret, s = 0, tmp = 0;

 /* for(n = 0; n < a->N; ++n){
    tmp = mu(a, n, argPhi, argAlp, argDel);
    s += a->y[iG(n, g)] * tmp - exp(a->c[iN(a->mC, n)] + 
         a->eps[iNG(a->mEps, n, g)] + tmp);
  } */
  
  s = 100;

  /* phi part */
  ret = s - pow(argPhi - a->thePhi[a->mThePhi], 2) / (2 * pow(a->sigPhi[a->mSigPhi], 2));

  /* alpha part */
  if(argAlp * argAlp > 1e-6){
    tmp = -pow(argAlp - a->theAlp[a->mTheAlp], 2) / (2 * pow(a->sigAlp[a->mSigAlp], 2)) -
                log(1 - a->piAlp[a->mPiAlp]);
  } else {
    tmp = log(a->piAlp[a->mPiAlp]);
  }

  ret = ret + tmp;

  /* delta part */
  if(argDel * argDel > 1e-6){
    tmp = -pow(argDel - a->theDel[a->mTheDel], 2) / (2 * pow(a->sigDel[a->mSigDel], 2)) -
                log(1 - a->piDel[a->mPiDel]);
  } else {
    tmp = log(a->piDel[a->mPiDel]);
  }

  return ret + tmp;
}

__global__ void samplePhiAlpDelJoint_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g = IDX, G = a->G;
  num_t oldPhi, newPhi, oldAlp, newAlp, oldDel, newDel;
  num_t dl, lp, lu;

  if(g < G){

    oldPhi = a->phi[iG(a->mPhi, g)];
    newPhi = rnormalDevice(a, g, oldPhi, a->tunePhi[g]);

    oldAlp = a->alp[iG(a->mAlp, g)];
    newAlp = alpProp(a, g);

    oldDel = a->del[iG(a->mDel, g)];
    newDel = delProp(a, g);

    dl =  lPhiAlpDelJoint(a, g, newPhi, newAlp, newDel) 
         -lPhiAlpDelJoint(a, g, oldPhi, oldAlp, oldDel);
    lp = 0 < dl ? 0 : dl;
    lu = log(runiformDevice(a, g, 0, 1));
    
    if(lu < lp){ /* accept */
      a->phi[iG(a->mPhi + 1, g)] = newPhi;
      a->alp[iG(a->mAlp + 1, g)] = newAlp;
      a->del[iG(a->mDel + 1, g)] = newDel;

      a->tunePhi[g] *= 1.1; 

      if(a->mPhi >= a->burnin){
        ++a->accPhi[g];
        ++a->accAlp[g];
        ++a->accDel[g];
      }
    } else { /* reject */
      a->phi[iG(a->mPhi + 1, g)] = oldPhi;
      a->alp[iG(a->mAlp + 1, g)] = oldAlp;
      a->del[iG(a->mDel + 1, g)] = oldDel;

      a->tunePhi[g] /= 1.1;
    }
  }
}

__global__ void samplePhiAlpDelJoint_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  ++a->mPhi;
  ++a->mAlp;
  ++a->mDel;
}

__host__ void samplePhiAlpDelJoint(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */
  float myTime;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  fprintf(cfg->log, "phiAlpDelJoint ");

  samplePhiAlpDelJoint_kernel1<<<G_GRID, G_BLOCK>>>(dev_a);
  samplePhiAlpDelJoint_kernel2<<<1, 1>>>(dev_a);

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  
  fprintf(cfg->time, "%0.3f ", myTime/MILLISECS); /* elapsed time */
}