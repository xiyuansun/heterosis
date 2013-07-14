#include <constants.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

num_t runiform(num_t lb, num_t ub){ 
  return (ub - lb) * (rand() / (num_t) RAND_MAX) + lb;
}