#ifndef CONFIG_H
#define CONFIG_H

#include <constants.h>

typedef struct {

  char *dataFile; /* space-delimited text file with G rows and N columns */
  char *groupFile; /* space-delimited text file with 1 row and N entries */
  char *probsFile; /* main output: diff expression / heterosis probabilities */
  char *ratesFile; /* optional output: acceptance rates of Metropolis steps */
  char *hyperFile; /* optional output: hyperparameters */
  char *parmsFile; /* optional output: example parameters */

  int probsFlag; /* indicate choice to output probs of diff expression / heterosis */
  int ratesFlag; /* indicate choice to output acceptance rates */
  int hyperFlag; /* indicate choice to output hyperparameters */
  int parmsFlag; /* indicate choice to output example parameters */
   
  int burnin; /* burn-in of chain. Defaults to M/2. */
  int joint; /* indicate joint sampling of phi_g, alpha_g, and delta_g */
  int seed; /* seed for random number generators */
  int heterosis; /* indicate whether the program will test for heterosis */

  int M; /* length of chain (not including initial values) */
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
  
  /* hyperparameters */
  
  float sigC;
  float d;
  float tau;
  float thePhi;
  float theAlp;
  float theDel;
  float sigPhi;
  float sigAlp;
  float sigDel;
  float piAlp;
  float piDel;
  
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