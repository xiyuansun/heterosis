#ifndef CONFIG_H
#define CONFIG_H

#include <constants.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {

  int chainNum; /* index of current chain */
  int m; /* current place in current chain */
  int chains; /* number of chains */

  char dataFile[BUF]; /* space-delimited text file with G rows and N columns */
  char groupFile[BUF]; /* space-delimited text file with 1 row and N entries */

  int ratesFlag; /* choice to output acceptance rates */
  int hyperFlag; /* choice to output hyperparameters */
  int parmsFlag; /* choice to output parameters */
  int timeFlag; /* output how much time it takes to sample each kind of parameter */
  int verbose;
  int diagnostics;
  
  /* curand states */
  
  curandState_t *states;
   
  int burnin; /* burn-in of chain. Defaults to M/2. */
  int joint; /* indicate joint sampling of phi_g, alpha_g, and delta_g */
  int seed; /* seed for random number generators */
  int heterosis; /* indicate whether the program will test for heterosis */

  int M; /* length of chain (not including initial values) */
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
  
  /* hyperparameters */
  
  num_t sigC;
  num_t d;
  num_t tau;
  num_t thePhi;
  num_t theAlp;
  num_t theDel;
  num_t sigPhi;
  num_t sigAlp;
  num_t sigDel;
  num_t piAlp;
  num_t piDel;
  
  /* indicate choices to hold each hyperparameter constant */
  
  int constSigC;
  int constD;
  int constTau;
  int constThePhi;
  int constTheAlp;
  int constTheDel;
  int constSigPhi;
  int constSigAlp;
  int constSigDel;
  int constPiAlp;
  int constPiDel;
  
  /* total time spent sampling each parameter */
  
  num_t timeC;
  num_t timeTau;
  num_t timePiAlp;
  num_t timePiDel;
  num_t timeD;
  num_t timeThePhi;
  num_t timeTheAlp;
  num_t timeTheDel;
  num_t timeSigC;
  num_t timeSigPhi;
  num_t timeSigAlp;
  num_t timeSigDel;
  num_t timeEta;
  num_t timeEps;
  num_t timePhi;
  num_t timeAlp;
  num_t timeDel;

} Config;

#endif /* CONFIG_H */