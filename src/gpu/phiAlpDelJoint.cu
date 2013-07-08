#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

num_t lPhiAlpDelJoint(Chain *a, int g, num_t argPhi, num_t argAlp, num_t argDel){ /* device */
 
  int n;
   N = a->N, G = a->G;
  
  num_t ret, s = 0, tmp = 0;

  for(n = 0; n < a->N; ++n){
    tmp = mu(a, n, argPhi, argAlp, argDel);
    s += a->y[iNG(n, g)] * tmp - exp(a->c[iMN(a->mC, n)] + 
         a->eps[iMNG(a->mEps, n, g)] + tmp);
  }

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

void samplePhiAlpDelJoint_kernel1(Chain *a){ /* kernel <<<G, 1>>> */
  int g;
  int G = a->G;
  
  num_t oldPhi, newPhi, oldAlp, newAlp, oldDel, newDel;
  num_t dl, lp, lu;

  for(g = 0; g < a->G; ++g){

    oldPhi = a->phi[iMG(a->mPhi, g)];
    newPhi = rnormal(oldPhi, a->tunePhi[g]);

    oldAlp = a->alp[iMG(a->mAlp, g)];
    newAlp = alpProp(a, g);

    oldDel = a->del[iMG(a->mDel, g)];
    newDel = delProp(a, g);

    dl = lPhiAlpDelJoint(a, g, newPhi, newAlp, newDel) 
       - lPhiAlpDelJoint(a, g, oldPhi, oldAlp, oldDel);
    lp = 0 < dl ? 0 : dl;
    lu = log(runiform(0, 1));
    
    if(lu < lp){ /* accept */
      a->phi[iMG(a->mPhi + 1, g)] = newPhi;
      a->alp[iMG(a->mAlp + 1, g)] = newAlp;
      a->del[iMG(a->mDel + 1, g)] = newDel;

      a->tunePhi[g] *= 1.1; 

      if(a->mPhi >= a->burnin){
        ++a->accPhi[g];
        ++a->accAlp[g];
        ++a->accDel[g];
      }
    } else { /* reject */
      a->phi[iMG(a->mPhi + 1, g)] = oldPhi;
      a->alp[iMG(a->mAlp + 1, g)] = oldAlp;
      a->del[iMG(a->mDel + 1, g)] = oldDel;

      a->tunePhi[g] /= 1.1;
    }
  }
}

void samplePhiAlpDelJoint_kernel2(Chain *a){ /* kernel <<<1, 1>>> */
  ++a->mPhi;
  ++a->mAlp;
  ++a->mDel;
}

void samplePhiAlpDelJoint(Chain *a){ /* host */
  samplePhiAlpDelJoint_kernel1(a);
  samplePhiAlpDelJoint_kernel2(a);
}