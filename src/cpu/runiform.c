#include <constants.h>
#include <math.h>
#include <stdlib.h>

num_t runiform(num_t lb, num_t ub){
  int max = 1e5;
  int i = rand() % max;
  num_t maxf = (num_t) max;
  num_t u = (num_t) i;
  
  u /= maxf;
  u = (ub - lb) * u + lb;
  return u;
}