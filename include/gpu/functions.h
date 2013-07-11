#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "Chain.h"
#include "Config.h"
#include "constants.h"
#include <cuda.h>
#include <cuda_runtime.h>

/* indexing macros */

#define iNG(m, n, g) ((m) * N * G + (n) * G + (g))
#define iN(m, n) ((m) * N + (n))
#define iG(n, g) ((n) * G + (g))

/* cuda call macros */

#define CUDA_CALL(x) {if((x) != cudaSuccess){ \
  printf("CUDA error at %s:%d\n",__FILE__,__LINE__); \
  printf("  %s\n", cudaGetErrorString(cudaGetLastError())); \
  exit(EXIT_FAILURE);}} 
  
#define CURAND_CALL(x) {if((x) != CURAND_STATUS_SUCCESS) { \
  printf("CURAND error at %s:%d\n",__FILE__,__LINE__); \
  exit(EXIT_FAILURE);}} 

#define FREE(x, onHost) {if(onHost){free(x);} else {CUDA_CALL(cudaFree(x));}}

/* function declarations */

__host__ void pi1(int*, int, const char*);
__host__ void pf1(num_t*, int, const char*);
__host__ void pi2(int*, int, int, const char*);
__host__ void pf2(num_t*, int, int, const char*);
__host__ void pstate(curandState*, int, int, const char*);
__host__ void pi3(int*, int, int, int, const char*);
__host__ void pf3(num_t*, int, int, int, const char*);

__host__ Config *config(int, char**);
__host__ void getopts(Config*, int, char**);
__host__ void printConfig(Config*);
__host__ void freeConfig(Config*);

__host__ int *mySampleInt(int, int);
__host__ int *readGrp(Config*);
__host__ count_t *readData(Config*);

__host__ void allocChainHost(Chain**, Config*);
__host__ void allocChainDevice(Chain**, Chain**, Config*);
  
__host__ int cmpfunc (const void*, const void*);
__global__ void curand_setup_kernel(Chain*, unsigned int);
                                  
__host__ void newChain(Chain** host_a, Chain **dev_a, Config*);
__global__ void newChain_kernel1(Chain*);
__global__ void newChain_kernel2(Chain*);
__host__ Chain *chainDeviceToHost(Chain*, Chain*, Config*);
__host__ void printChain(Chain*, Chain*, Config*);
__host__ void freeChain(Chain*, Config*, int);

__host__ num_t runiform(num_t, num_t);
__host__ num_t rnormal(num_t, num_t);
__host__ num_t rgamma(num_t, num_t, num_t);
__host__ num_t rbeta(num_t, num_t);

__global__ void lC_kernel1(Chain*, int);
__global__ void lC_kernel2(Chain*, int, int);
__host__ void lC(Chain*, Chain*, Config*, int, int);
__global__ void sampleC_kernel1(Chain*);
__global__ void sampleC_kernel2(Chain*);
__global__ void sampleC_kernel3(Chain*);
__host__ void sampleC(Chain*, Chain*, Config*);

__device__ num_t lEps(Chain*, int, int, num_t);
__global__ void sampleEps_kernel1(Chain*);
__global__ void sampleEps_kernel2(Chain*);
__host__ void sampleEps(Chain*, Chain*, Config*);

__global__ void lD_kernel1(Chain*, int);
__global__ void lD_kernel2(Chain*);
__global__ void lD_kernel3(Chain*, int);
__host__ void lD(Chain*, Chain*, Config *cfg, int);
__global__ void sampleD_kernel1(Chain*);
__global__ void sampleD_kernel2(Chain*);
__host__ void sampleD(Chain*, Chain*, Config*);

num_t lPhi(Chain*, int, num_t);
void samplePhi_kernel1(Chain*);
void samplePhi_kernel2(Chain*);
void samplePhi(Chain*, Chain*, Config*);

num_t alpProp(Chain*, int);
num_t lAlp(Chain*, int, num_t);
void sampleAlp_kernel1(Chain*);
void sampleAlp_kernel2(Chain*);
void sampleAlp(Chain*, Chain*, Config*);

num_t delProp(Chain*, int);
num_t lDel(Chain*, int, num_t);
void sampleDel_kernel1(Chain*);
void sampleDel_kernel2(Chain*);
void sampleDel(Chain*, Chain*, Config*);

num_t lPhiAlpDelJoint(Chain*, int, num_t, num_t, num_t);
void samplePhiAlpDel_kernel1Joint(Chain*);
void samplePhiAlpDel_kernel2Joint(Chain*);
void samplePhiAlpDelJoint(Chain*, Chain*, Config*);

void samplePhiAlpDel(Chain*, Chain*, Config*);

__global__ void sampleSigC_kernel(Chain*);
__host__ void sampleSigC(Chain*, Chain*, Config*);

__global__ void sampleEta_kernel1(Chain*); 
__global__ void sampleEta_kernel2(Chain*);
__global__ void sampleEta_kernel3(Chain*);
__host__ void sampleEta(Chain*, Chain*, Config*);

__global__ void sampleTau_kernel1(Chain*);
__global__ void sampleTau_kernel2(Chain*);
__host__ void sampleTau(Chain*, Chain*, Config*);

__global__ void sampleThePhi_kernel1(Chain*);
__global__ void sampleThePhi_kernel2(Chain*);
__host__ void sampleThePhi(Chain*, Chain*, Config*);

__global__ void sampleTheAlp_kernel1(Chain*);
__global__ void sampleTheAlp_kernel2(Chain*);
__host__ void sampleTheAlp(Chain*, Chain*, Config*);

__global__ void sampleTheDel_kernel1(Chain*);
__global__ void sampleTheDel_kernel2(Chain*);
__host__ void sampleTheDel(Chain*, Chain*, Config*);

__global__ void sampleSigPhi_kernel1(Chain*);
__global__ void sampleSigPhi_kernel2(Chain*);
__host__ void sampleSigPhi(Chain*, Chain*, Config*);

__global__ void sampleSigAlp_kernel1(Chain*);
__global__ void sampleSigAlp_kernel2(Chain*);
__host__ void sampleSigAlp(Chain*, Chain*, Config*);

__global__ void sampleSigDel_kernel1(Chain*);
__global__ void sampleSigDel_kernel2(Chain*);
__host__ void sampleSigDel(Chain*, Chain*, Config*);

__global__ void samplePiAlp_kernel1(Chain*);
__global__ void samplePiAlp_kernel2(Chain*);
__host__ void samplePiAlp(Chain*, Chain*, Config*);

__global__ void samplePiDel_kernel1(Chain*);
__global__ void samplePiDel_kernel2(Chain*);
__host__ void samplePiDel(Chain*, Chain*, Config*);

__host__ void runChain(Chain*, Chain*, Config*);
__host__ void oneChain(int, char**);

__host__ void printProbs(Chain*, Config*);
__host__ void printRates(Chain*, Config*);
__host__ void printHyper(Chain*, Config*);
__host__ void printParms_oneFile(Chain*, Config*, int);
__host__ void printParms(Chain*, Config*);
__host__ void summarizeChain(Chain*, Chain*, Config*);

/* definitions of inline functions */

inline __device__ num_t mu(Chain *a, int n, num_t phi, num_t alp, num_t del){
  if(a->grp[n] == 1){
    return phi - alp;
  } else if(a->grp[n] == 2){
    return phi + del;
  } else {
    return phi + alp;
  }
}

inline __device__ num_t runiformDevice(Chain *a, int g, num_t lb, num_t ub){
  num_t u = curand_uniform(&(a->states[g]));
  return (ub - lb) * u + lb;
}

inline __device__ num_t rnormalDevice(Chain *a, int g, num_t m, num_t s){

  num_t u1 = runiformDevice(a, g, 0, 1);
  num_t u2 = runiformDevice(a, g, 0, 1);
  
  return sqrt(-2 * log(u1)) * sin(2 * M_PI * u2) * s + m;
}

inline __device__ num_t rgammaDevice(Chain *a, int g, num_t shape, num_t rate, num_t lb){
   
  num_t A, c, d, r, u, v, w, x, z, eps, eps0, 
        haznaz, lam, ret, tmp1, tmp2;

  if(shape <= 0){
    return(-1);
  }
  
  if(rate <= 0){
    return(-1);
  }

  if((shape >= 1) && (shape - 1 < lb * rate)){ /* Chung (1998) */

    c = lb * rate;
    eps0 = (c - shape + sqrt(pow(c - shape, 2) + 4 * c))/(2 * c);

    if(eps0 > 1){
      eps = 0.75;
    } else {
      eps = eps0;
    }

    A = pow((1 - eps)/(shape - 1), shape - 1) * exp(shape - 1);

    while(1){
      x = - (1/eps) * log(runiformDevice(a, g, 0, 1)) + lb * rate;
      u = runiformDevice(a, g, 0, 1);

      if(u < A * pow(x, shape - 1) * exp((eps - 1) * x))
        return(x / rate);
    }

  } else if(shape >= 1){ /* Marsaglia and Tsang (2000) */

    d = shape - 1/3;
    c = 1 / sqrt(9 * d);

    while(1){
      v = -1;
      while(v <= 0){
        x = rnormalDevice(a, g, 0, 1);
        v = pow(1 + c*x, 3);
      }

      ret = d * v / rate;

      if(ret > lb){
        u = runiformDevice(a, g, 0, 1);

        if(u < 1 - 0.0331 * pow(x, 4))
          return(ret);

        if(log(u) < 0.5 * pow(x, 2) + d * (1 - v + log(v)))
          return(ret);
      }
    }
  } else if (0.135 <= shape && shape < 1){ /* Kundu and Gupta (2006) */

    while(1){      

      u = runiformDevice(a, g, 0, 1);
      x = -2 * log(1 - pow(u, 1 / shape));
      ret = x / rate;

      if(ret > lb){
        v = runiformDevice(a, g, 0, 1);

        tmp1 = exp(-x/2);
        tmp2 = pow(x, shape - 1)* tmp1 * pow(((float) 2), ((float) (1 - shape))) * pow(((float) (1 - tmp1)), ((float) (1 - shape)));

        if(v < tmp2)
          return(ret);
      }
    }
  } else{ /* Martin and Liu (2013) */
   
    while(1){  
      lam = 1/shape - 1;
      w = shape / (exp(1 - shape));
      r = 1 / (1 + w); 
      u = runiformDevice(a, g, 0, 1);

      if(u <= r){
        z = - log(u / r);
      } else {
        z = log(runiformDevice(a, g, 0, 1)) / lam;
      }
      
      ret = exp(-z / shape) / rate;

      if(ret > lb){
        if(z >= 0){
          haznaz = exp(-exp(-z / shape));
        } else{
          haznaz = 1 / (w * lam) * exp((lam - 1) * z - exp(z / shape));
        }

        if(haznaz > runiformDevice(a, g, 0, 1))
          return(ret);
      }
    }
  }
}

inline __device__ num_t rbetaDevice(Chain *chain, int g, num_t a, num_t b){

  num_t x = rgammaDevice(chain, g, a, 1, 0);
  num_t y = rgammaDevice(chain, g, b, 1, 0);
  
  return x / (x + y);
}

#endif /* FUNCTIONS_H */
