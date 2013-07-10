#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdlib.h>

__host__ num_t rbeta(num_t a, num_t b){

  num_t x = rgamma(a, 1, 0);
  num_t y = rgamma(b, 1, 0);
  
  return x / (x + y);
}

inline __device__ num_t rbetaDevice(Chain *chain, int g, num_t a, num_t b){

  num_t x = rgammaDevice(chain, g, a, 1, 0);
  num_t y = rgammaDevice(chain, g, b, 1, 0);
  
  return x / (x + y);
}