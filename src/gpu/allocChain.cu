#include <Chain.h>
#include <Config.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>


__global__ void allocChain_kernel(Chain *a, num_t *y, num_t *yMeanG, int *grp, num_t *c, num_t *sigC,
                                  num_t *eps, num_t *eta, num_t *d, num_t *tau, num_t *phi, num_t *thePhi,
                                  num_t *sigPhi, num_t *alp, num_t *theAlp, num_t *sigAlp, 
                                  num_t *piAlp, num_t *del, num_t *theDel, num_t *sigDel, 
                                  num_t *piDel, num_t *tmp1, num_t *tmp2, num_t *Old, num_t *Nw, 
                                  num_t *lOld, num_t *lNew, num_t *tuneC, num_t *tunePhi,
                                  num_t *tuneEps, int *accC, int *accPhi, int *accAlp, 
                                  int *accDel, int *accEps){

  a->y = y;
  a->yMeanG = yMeanG;
  a->grp = grp;
    
  a->c = c;
  a->sigC = sigC;
  a->eps = eps;
  a->eta = eta;
  a->d = d;
  a->tau = tau;
  a->phi = phi;
  a->thePhi = thePhi;
  a->sigPhi = sigPhi;
  a->alp = alp;
  a->theAlp = theAlp;
  a->sigAlp = sigAlp;
  a->piAlp = piAlp;
  a->del = del;
  a->theDel = theDel;
  a->sigDel = sigDel;
  a->piDel = piDel;
    
  a->tmp1 = tmp1;
  a->tmp2 = tmp2;
  a->Old = Old;
  a->New = Nw;
  a->lOld = lOld;
  a->lNew = lNew;
    
  a->tuneC = tuneC;
  a->tunePhi = tunePhi;
  a->tuneEps = tuneEps;
    
  a->accC = accC;
  a->accPhi = accPhi;
  a->accAlp = accAlp;
  a->accDel = accDel;
  a->accEps = accEps;
}


__host__ Chain *allocChain(Config *cfg, int onHost){

  Chain *a;
  
  int *grp, *accC, *accPhi, *accAlp, *accDel, *accEps;
  num_t *y, *yMeanG, *c, *sigC, *eps, *eta, *d, *tau, *phi, *thePhi *sigPhi, *alp, 
        *theAlp, *sigAlp, *piAlp, *del, *theDel, *sigDel, *piDel;
  num_t *tmp1, *tmp2, *Old, *Nw, *lOld, *lNew;
  num_t *tuneC, *tunePhi, *tuneEps;
  
  /* data */  
      
  ALLOC(y, (count_t*), cfg->N * cfg->G * sizeof(count_t), onHost);
  ALLOC(yMeanG, (num_t*), cfg->N * sizeof(num_t), onHost);
  ALLOC(grp, (int*), cfg->N * sizeof(int), onHost);

  /* parameters */
  
  ALLOC(c, (num_t*), (cfg->M + 1) * cfg->N * sizeof(num_t), onHost);
  ALLOC(sigC, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(eps, (num_t*), (cfg->M + 1) * cfg->N * cfg->G * sizeof(num_t), onHost);
  ALLOC(eta, (num_t*), (cfg->M + 1) * cfg->G * sizeof(num_t), onHost);
  ALLOC(d, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(tau, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(phi, (num_t*), (cfg->M + 1) * cfg->G * sizeof(num_t), onHost);
  ALLOC(thePhi, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(sigPhi, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(alp, (num_t*), (cfg->M + 1) * cfg->G * sizeof(num_t), onHost);
  ALLOC(theAlp, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(sigAlp, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(piAlp, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(del, (num_t*), (cfg->M + 1) * cfg->G * sizeof(num_t), onHost);
  ALLOC(theDel, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(sigDel, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  ALLOC(piDel, (num_t*), (cfg->M + 1) * sizeof(num_t), onHost);
  
  /* temporary and return values */
    
  ALLOC(tmp1, (num_t*), cfg->G * sizeof(num_t), onHost);
  ALLOC(tmp2, (num_t*), cfg->G * sizeof(num_t), onHost);
  
  ALLOC(Old, (num_t*), cfg->N * sizeof(num_t), onHost);
  ALLOC(New, (num_t*), cfg->N * sizeof(num_t), onHost);
  ALLOC(lOld, (num_t*), cfg->N * sizeof(num_t), onHost);
  ALLOC(lNew, (num_t*), cfg->N * sizeof(num_t), onHost);

  /* tuning parameters for Metropolis steps */
  
  ALLOC(tuneC, (num_t*), cfg->N * sizeof(num_t), onHost);
  ALLOC(tunePhi, (num_t*), cfg->G * sizeof(num_t), onHost);
  ALLOC(tuneEps, (num_t*), cfg->N * cfg->G * sizeof(num_t), onHost);

  /* number of acceptances for Metropolis steps */
  
  ALLOC(accC, (int*), cfg->N * sizeof(int), onHost);
  ALLOC(accPhi, (int*), cfg->G * sizeof(int), onHost);
  ALLOC(accAlp, (int*), cfg->G * sizeof(int), onHost);
  ALLOC(accDel, (int*), cfg->G * sizeof(int), onHost);
  ALLOC(accEps, (int*), cfg->N * cfg->G * sizeof(int), onHost);
   
  if(onHost){
    a = (Chain*) malloc(sizeof(Chain));
    
    a->y = y;
    a->yMeanG = yMeanG;
    a->grp = grp;
    
    a->c = c;
    a->sigC = sigC;
    a->eps = eps;
    a->eta = eta;
    a->d = d;
    a->tau = tau;
    a->phi = phi;
    a->thePhi = thePhi;
    a->sigPhi = sigPhi;
    a->alp = alp;
    a->theAlp = theAlp;
    a->sigAlp = sigAlp;
    a->piAlp = piAlp;
    a->del = del;
    a->theDel = theDel;
    a->sigDel = sigDel;
    a->piDel = piDel;
    
    a->tmp1 = tmp1;
    a->tmp2 = tmp2;
    a->Old = Old;
    a->New = Nw;
    a->lOld = lOld;
    a->lNew = lNew;
    
    a->tuneC = tuneC;
    a->tunePhi = tunePhi;
    a->tuneEps = tuneEps;
    
    a->accC = accC;
    a->accPhi = accPhi;
    a->accAlp = accAlp;
    a->accDel = accDel;
    a->accEps = accEps;
    
  } else {
    CUDA_CALL(cudaMalloc((void **) &a, sizeof(Chain)));
    allocChain_kernel<<<1, 1>>>(a, y, yMeanG, grp, c, sigC, eps, eta, d, tau, phi, thePhi,
                                sigPhi, alp, theAlp, sigAlp, piAlp, del, theDel, sigDel, 
                                piDel, tmp1, tmp2, Old, Nw, lOld, lNew, tuneC, tunePhi,
                                tuneEps, accC, accPhi, accAlp, accDel, accEps);
  }
    
  return a;
}