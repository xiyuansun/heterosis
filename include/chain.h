#ifndef CHAIN_H
#define

#define typedef struct {

  int **y; /* data */
  int *yMeanG; /* gene-wise mean expression levels */
  int *grp; /* assignment of libraries (samples) to treatment groups */ 

  int M; /* length of chain */
  int N; /* number of libraries (samples) */
  int G; /* number of genes */
  
  /* initialization constants */
  
  float sigC0;
  float d0;
  float aTau;
  float aAlp;
  float aDel;
  float bTau;
  float bAlp;
  float bDel;
  float gamPhi;
  float gamAlp;
  float gamDel;
  float sigPhi0;
  float sigAlp0;
  float sigDel0;

  /* samples of parameters and hyperparameters */
  
  float **c;
    float *sigC;
   
  float ***eps;
    float **eta;
      float *d;
      float *tau;
      
  float **phi;
    float *thePhi;
    float *sigPhi;
    
  float **alp;
    float *theAlp;
    float *sigAlp;
    float *piAlp;
    
  float **del;
    float *theDel;
    float *sigDel;
    float *piDel;
    
  /* temporary and return values of kernels */

  float shape;
  float rate;
  
  float *tmp1;
  float *tmp2;
  
  float **Old;
  float **New;
  
  float **lOld;
  float **lNew;
  
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
  
  float tunC;
  float tunEps;
  float tunD;
  float tunPhi;
  
  /* number of acceptances for metropolis steps */

  int accD;
  int *accC;
  int **accEps;
  int *accPhi;
  int *accAlp;
  int *accDel;

} Chain;

#endif /* CHAIN_H */