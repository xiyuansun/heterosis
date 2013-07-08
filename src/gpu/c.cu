#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void lC_kernel1(Chain *a, int n){ /* kernel <<<G, 1>>> */
  int g;
  int M = a->M, N = a->N, G = a->G;
  
  for(g = 0; g < G; ++g)
    a->tmp1[g] = exp(a->eps[iMNG(a->mEps, n, g)] + mu(a, n, a->phi[iMG(a->mPhi, g)], 
                    a->alp[iMG(a->mAlp, g)], a->del[iMG(a->mDel, g)]));
}

void lC_kernel2(Chain *a, int n){ /* parallel pairwise sum */
  int g;
  int M = a->M, N = a->N, G = a->G;;
  a->tmp2[0] = 0;
  
  for(g = 0; g < G; ++g)
    a->tmp2[0] += a->tmp1[g];
}

void lC_kernel3(Chain *a, int n, int newArg){ /* kernel <<<1, 1>>> */
  num_t arg, ret;
  int M = a->M, N = a->N, G = a->G;

  if(newArg){
    arg = a->New[n];
  } else {
    arg = a->Old[n];
  }

  ret = arg * G * a->yMeanG[n] - exp(arg) * a->tmp2[0] - (arg*arg) / 
        (2 * a->sigC[a->mSigC] * a->sigC[a->mSigC]);

  if(newArg){
    a->lNew[n] = ret;
  } else {
    a->lOld[n] = ret;
  }
}

void lC(Chain *a, int n, int newArg){ /* host */
  lC_kernel1(a, n);
  lC_kernel2(a, n);
  lC_kernel3(a, n, newArg);
}

void sampleC_kernel1(Chain *a){ /* kernel <<<1, N>>> */
  int n;
  int M = a->M, N = a->N, G = a->G;
  
  for(n = 0; n < N; ++n){
    a->Old[n] = a->c[iMN(a->mC, n)];
    a->New[n] = rnormal(a->Old[n], a->tuneC[n]);
  }
}

void sampleC_kernel2(Chain *a){ /* kernel <<<1, N>>> */
  int n;
  int M = a->M, N = a->N, G = a->G;
  num_t dl, lp, lu;

  for(n = 0; n < N; ++n){ 

    dl = a->lNew[n] - a->lOld[n];
    lp = 0 < dl ? 0 : dl;
    lu = log(runiform(0, 1));
      
    if(lu < lp){ /* accept */
      a->c[iMN(a->mC + 1, n)] = a->New[n];
      a->tuneC[n] *= 1.1; /* Increase the proposal variance to avoid  
                                       getting stuck in a mode */
                                       
      if(a->mC >= a->burnin)                                 
        ++a->accC[n];
    } else { /* reject */
      a->c[iMN(a->mC + 1, n)] = a->Old[n];
      a->tuneC[n] /= 1.1; /* If you're rejecting too often, decrease the  
                                       proposal variance to sample closer to 
                                       the last accepted value. */
    }
  }
}

void sampleC_kernel3(Chain *a){ /* kernel <<<1, 1>>> */
  ++a->mC;
}

void sampleC(Chain *a){ /* host */
  int n;
  int M = a->M, N = a->N, G = a->G;
  sampleC_kernel1(a);

  for(n = 0; n < N; ++n){ 
    lC(a, n, 1);
    lC(a, n, 0);
  }

  sampleC_kernel2(a);
  sampleC_kernel3(a);
}