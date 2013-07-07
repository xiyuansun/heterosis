#ifndef SUMMARY_H
#define SUMMARY_H

#include <constants.h>

typedef struct {

  /* hyperparameters */

  num_t *sigC;
  num_t *d;
  num_t *tau;
  num_t *thePhi;
  num_t *theAlp;
  num_t *theDel;
  num_t *sigPhi;
  num_t *sigAlp;
  num_t *sigDel;
  num_t *piAlp;
  num_t *piDel;
  
  /* acceptance rates of metropolis steps */
  
  num_t accD;
  num_t *accC;
  num_t *accEps;
  num_t *accPhi;
  num_t *accAlp;
  num_t *accDel;
  
  /* probabilities of differential expression and heterosis */
  
  num_t *prob_de;
  num_t *prob_hph;
  num_t *prob_lph;
  num_t *prob_mph;
  
  /* genes and libraries from which to take example parameters */
  
  int nlibs;
  int ngenes;
  
  int *libs;
  int *genes;
  
  /* parameters */
  
  num_t **c;
  num_t ***eps;
  num_t **eta;
  num_t **phi;
  num_t **alp;
  num_t **del;
  
} Summary;

#endif /* SUMMARY_H */