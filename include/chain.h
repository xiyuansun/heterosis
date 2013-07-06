#ifndef CHAIN_H
#define CHAIN_H

#include <numericTypes.h>

typedef struct {

  count_t **y; /* data */
  num_t *yMeanG; /* gene-wise mean expression levels */
  int *grp; /* assignment of libraries (samples) to treatment groups */ 

  int M; /* length of chain */
  int N; /* number of libraries (samples) */
  int G; /* number of genes */
  
  /* initialization constants */
  
  num_t sigC0;
  num_t d0;
  num_t aTau;
  num_t aAlp;
  num_t aDel;
  num_t bTau;
  num_t bAlp;
  num_t bDel;
  num_t gamPhi;
  num_t gamAlp;
  num_t gamDel;
  num_t sigPhi0;
  num_t sigAlp0;
  num_t sigDel0;

  /* samples of parameters and hyperparameters */
  
  num_t **c;
    num_t *sigC;
   
  num_t ***eps;
    num_t **eta;
      num_t *d;
      num_t *tau;
      
  num_t **phi;
    num_t *thePhi;
    num_t *sigPhi;
    
  num_t **alp;
    num_t *theAlp;
    num_t *sigAlp;
    num_t *piAlp;
    
  num_t **del;
    num_t *theDel;
    num_t *sigDel;
    num_t *piDel;
    
  /* temporary and return values of kernels */

  num_t shape;
  num_t rate;
  
  num_t *tmp1;
  num_t *tmp2;
  
  num_t **Old;
  num_t **New;
  
  num_t **lOld;
  num_t **lNew;
  
  /* current place in the chain of each parameter */
  
  int mC;
    int mSigC;
    
  int mEps;
    int mEta;
      int mD;
      int mTau;
      
  int mPhi;
    int mThePhi;
    int mSigPhi;
    
  int mAlp;
    int mTheAlp;
    int mSigAlp;
    int mPiAlp;
    
  int mDel;
    int mTheDel;
    int mSigDel;
    int mPiDel;
    
  /* tuning parameters for metropolis steps */
  
  num_t *tuneC;
  num_t **tuneEps;
  num_t tuneD;
  num_t *tunePhi;
  
  /* number of acceptances for metropolis steps */

  int accD;
  int *accC;
  int **accEps;
  int *accPhi;
  int *accAlp;
  int *accDel;

} Chain;

#endif /* CHAIN_H */