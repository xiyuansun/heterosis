#ifndef CONFIG_H
#define CONFIG_H

#include "constants.h"

typedef struct {

  int chainNum;

  char dataFile[BUF]; /* space-delimited text file with G rows and N columns */
  char groupFile[BUF]; /* space-delimited text file with 1 row and N entries */
  char logFile[BUF];
  
  int ratesFlag; /* choice to output acceptance rates */
  int hyperFlag; /* choice to output hyperparameters */
  int parmsFlag; /* choice to output parameters */
   
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

} Config;

#endif /* CONFIG_H */