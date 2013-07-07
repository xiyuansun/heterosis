#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void lC_kernel1(Chain *a, int n){ /* kernel <<<G, 1>>> */
  int g;
  
  for(g = 0; g < a->G; ++g)
    a->tmp1[g] = exp(a->eps[a->mEps][n][g] + mu(a, n, a->phi[a->mPhi][g], 
                    a->alp[a->mAlp][g], a->del[a->mDel][g]));
}

void lC_kernel2(Chain *a, int n){ /* parallel pairwise sum */
  int g;
  a->tmp2[0] = 0;
  
  for(g = 0; g < a->G; ++g)
    a->tmp2[0] += a->tmp1[g];
}

void lC_kernel3(Chain *a, int n, int newArg){ /* kernel <<<1, 1>>> */
  num_t arg, ret;

  if(newArg){
    arg = a->New[0][n];
  } else {
    arg = a->Old[0][n];
  }

  ret = arg * a->G * a->yMeanG[n] - exp(arg) * a->tmp2[0] - (arg*arg) / 
        (2 * a->sigC[a->mSigC] * a->sigC[a->mSigC]);

  if(newArg){
    a->lNew[0][n] = ret;
  } else {
    a->lOld[0][n] = ret;
  }
}

void lC(Chain *a, int n, int newArg){ /* host */
  lC_kernel1(a, n);
  lC_kernel2(a, n);
  lC_kernel3(a, n, newArg);
}

void sampleC_kernel1(Chain *a){ /* kernel <<<1, N>>> */
  int n;
  
  for(n = 0; n < a->N; ++n){
    a->Old[0][n] = a->c[a->mC][n];
    a->New[0][n] = normalHost(a->Old[0][n], a->tuneC[n]);
  }
}

void sampleC_kernel2(Chain *a){ /* kernel <<<1, N>>> */
  int n;
  num_t dl, lp, lu;

  for(n = 0; n < a->N; ++n){ 

    dl = a->lNew[0][n] - a->lOld[0][n];
    lp = 0 < dl ? 0 : dl;
    lu = log(uniformHost(0, 1));
      
    if(lu < lp){ /* accept */
      a->c[a->mC + 1][n] = a->New[0][n];
      a->tuneC[n] = a->tuneC[n] * 1.1; /* Increase the proposal variance to avoid  
                                       getting stuck in a mode */
                                       
      if(a->mC >= a->burnin)                                 
        ++a->accC[n];
    } else { /* reject */
      a->c[a->mC + 1][n] = a->Old[0][n];
      a->tuneC[n] = a->tuneC[n] / 1.1; /* If you're rejecting too often, decrease the  
                                       proposal variance to sample closer to 
                                       the last accepted value. */
    }
  }
}

void sampleC_kernel3(Chain *a){ /* kernel <<<1, 1>>> */
  a->mC = a->mC + 1;
}

void sampleC(Chain *a){ /* host */
  int n;
  sampleC_kernel1(a);

  for(n = 0; n < a->N; ++n){ 
    lC(a, n, 1);
    lC(a, n, 0);
  }

  sampleC_kernel2(a);
  sampleC_kernel3(a);
}