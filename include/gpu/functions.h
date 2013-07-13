#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <cuda.h>
#include <curand_kernel.h>

#define iG(n, g) ((n) * G + (g))

#define CUDA_CALL(x) {}
  
#define FREE(x, onHost) {if(onhost){free(x);} else {CUDA_CALL(cudaFree(x))}}

__host__ void pi1(int*, int, const char*);
__host__ void pf1(num_t*, int, const char*);
__host__ void pstate(curandState_t*, int, const char*);
__host__ void pi2(int*, int, int, const char*);
__host__ void pf2(num_t*, int, int, const char*);
__host__ void pi3(int*, int, int, int, const char*);
__host__ void pf3(num_t*, int, int, int, const char*);

__host__ Config *config(int, char**);
__host__ void getopts(Config*, int, char**);
__host__ void printConfig(Config*);

__host__ int *readGrp(Config*);
__host__ count_t *readData(Config*);

__host__ Chain *allocChainHost(Config*);
__host__ void allocChainDevice(Chain**, Chain**, Config*);
__host__ Chain *chainDeviceToHost(Chain*, Chain*, Config*);

__host__ void cmpfunc(const void*, const void*);
__global__ void curand_setup_kernel(Chain*, int*);

__global__ void newChain_kernel1(Chain*);
__global__ void newChain_kernel2(Chain*);
__host__ void newChain(Chain**, Chain**, Config*);
__host__ void printChain(Chain*, Chain*, Config*);
__host__ void freeChain(Chain*, Config*);

__host__ num_t runiform(num_t, num_t);
__host__ num_t rnormal(num_t, num_t);
__host__ num_t rgamma(num_t, num_t, num_t);
__host__ num_t rbeta(num_t, num_t);

__global__ void lC_kernel1(Chain*, int);
__global__ void lC_kernel2(Chain*, int, int);
__host__ void lC(Chain*, Chain*, int, int);
__global__ void sampleC_kernel1(Chain*);
__global__ void sampleC_kernel2(Chain*);
__host__ void sampleC(Chain*, Chain*, Config*);

__device__ num_t lEps(Chain*, int, int, num_t);
__global__ void sampleEps_kernel(Chain*);
__host__ void sampleEps(Chain*, Chain*, Config*);

__global__ void lD_kernel1(Chain*);
__global__ void lD_kernel2(Chain*, int);
__host__ void lD(Chain*, int);
__global__ void sampleD_kernel1(Chain*);
__global__ void sampleD_kernel2(Chain*);
__host__ void sampleD(Chain*, Chain*, Config*);

__device__ num_t lPhi(Chain*, int, num_t);
__global__ void samplePhi_kernel(Chain*);
__host__ void samplePhi(Chain*, Chain*, Config*);

__device__ num_t lAlp(Chain*, int, num_t);
__global__ void sampleAlp_kernel(Chain*);
__host__ void sampleAlp(Chain*, Chain*, Config*);


__device__ num_t lDel(Chain*, int, num_t);
__global__ void sampleDel_kernel(Chain*);
__host__ void sampleDel(Chain*, Chain*, Config*);

__device__ num_t lPhiAlpDelJoint(Chain*, int, num_t, num_t, num_t);
__global__ void samplePhiAlpDelJoint_kernel(Chain*);
__host__ void samplePhiAlpDelJoint(Chain*, Chain*, Config*);

__host__ void samplePhiAlpDel(Chain*, Chain*, Config*);

__global__ void sampleSigC_kernel(Chain*);
__host__ void sampleSigC(Chain*, Chain*, Config*);

__global__ void sampleEta_kernel1(Chain*);
__global__ void sampleEta_kernel2(Chain*);
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
__host__ void oneChain(Config*);
__host__ void chains(int, char**);

__host__ void printHeaders(Chain*, Chain*, Config*);
__host__ void interimResults(Chain*, Chain*, Config*);
__host__ void summarizeChain(Chain*, Chain*, Config*);

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
   
  num_t lA, c, d, r, u, lu, v, w, x, z, eps, eps0, 
        haznaz, lam, ret, tmp1, tmp2, n, nmax = 100;
  
  if(shape <= 0 || rate <= 0)
    return(-1);

  if((shape >= 1) && (shape - 1 < lb * rate)){ /* Chung (1998) */

    c = lb * rate;
    eps0 = (c - shape + sqrt(pow(c - shape, 2) + 4 * c))/(2 * c);

    if(eps0 > 1){
      eps = 0.75;
    } else {
      eps = eps0;
    }

    tmp1 = shape - 1;
    lA = tmp1 * (log(1 - eps) - log(tmp1)) + tmp1;

    for(n = 0; n < nmax; ++n){
      x = - (1/eps) * log(runiformDevice(a, g, 0, 1)) + lb * rate;
      lu = log(runiformDevice(a, g, 0, 1));

      if(lu < lA + tmp1 * log(x) + (eps - 1) * x)
        return (x / rate);
    }
    return (x / rate);

  } else if(shape >= 1){ /* Marsaglia and Tsang (2000) */

    d = shape - 1/3;
    c = 1 / sqrt(9 * d);

    for(n = 0; n < nmax; ++n){
      v = -1;
      while(v <= 0){
        x = rnormal(0, 1);
        v = pow(1 + c*x, 3);
      }

      ret = d * v / rate;

      if(ret > lb){
        u = runiformDevice(a, g, 0, 1);

        if(u < 1 - 0.0331 * pow(x, 4))
          return(ret);

        if(log(u) < 0.5 * pow(x, 2) + d * (1 - v + log(v)))
          return ret;
      }
    }
    return ret;
    
  } else if (0.135 <= shape && shape < 1){ /* Kundu and Gupta (2006) */

    for(n = 0; n < nmax; ++n){      

      u = runiformDevice(a, g, 0, 1);
      x = -2 * log(1 - pow(u, 1 / shape));
      ret = x / rate;

      if(ret > lb){
        v = runiformDevice(a, g, 0, 1);

        tmp1 = exp(-x/2);
        tmp2 = pow(x, shape - 1)* tmp1 * pow(2, 1 - shape) * pow(1 - tmp1, 1 - shape);

        if(v < tmp2)
          return ret;
      }
    }
    return ret;
    
  } else{ /* Martin and Liu (2013) */
   
    for(n = 0; n < nmax; ++n){  
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
          return ret;
      }
    }
    
    return ret;
  }
}

inline __device__ num_t rbetaDevice(Chain *chain, int g, num_t a, num_t b){

  num_t x = rgammaDevice(chain, g, a, 1, 0);
  num_t y = rgammaDevice(chain, g, b, 1, 0);
  
  return x / (x + y);
}

inline __device__ num_t alpProp(Chain *a, int g){ /* device */

  num_t gam = a->gamAlp;
  num_t sig = a->sigAlp;

  num_t gprec = 1/(gam * gam);
  num_t sprec = 1/(sig * sig);

  num_t avg = (a->alp[g] * sprec) / (gprec + sprec);
  num_t s = gam * gam + sig * sig;
  num_t u = runiformDevice(a, g, 0, 1); 
  num_t nw;
  
  if(u < a->piAlp){
    nw = 0;
  } else {
    nw = rnormalDevice(a, g, avg, s);
  }

  return nw;
}

inline __device__ num_t delProp(Chain *a, int g){ /* device */   

  num_t gam = a->gamDel;
  num_t sig = a->sigDel;

  num_t gprec = 1/(gam * gam);
  num_t sprec = 1/(sig * sig);

  num_t avg = (a->del[g] * sprec) / (gprec + sprec);
  num_t s = gam * gam + sig * sig;
  num_t u = runiformDevice(a, g, 0, 1);
  num_t nw;

  if(u < a->piDel){
    nw = 0;
  } else {
    nw = rnormalDevice(a, g, avg, s);
  }

  return nw;
}

#endif /* FUNCTIONS_H */