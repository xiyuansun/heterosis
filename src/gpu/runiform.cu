#include <Chain.h>
#include <constants.h>
#include <curand_kernel.h>
#include <functions.h>
#include <math.h>
#include <stdlib.h>

__host__ num_t runiform(num_t lb, num_t ub){
  num_t u = ((num_t) rand()) / ((num_t) RAND_MAX);
  return (ub - lb) * u + lb;
}

__device__ num_t runiformDevice(Chain *a, int g, num_t lb, num_t ub){
  num_t u = curand_uniform(&(a->states[g]));
  return (ub - lb) * u + lb;
}