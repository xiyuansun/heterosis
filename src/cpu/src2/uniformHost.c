#include <constants.h>
#include <math.h>
#include <stdlib.h>

num_t uniformHost(num_t lb, num_t ub){
  num_t u = ((num_t) rand()) / ((num_t) RAND_MAX);
  return (ub - lb) * u + lb;
}