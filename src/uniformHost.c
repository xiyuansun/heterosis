#include <math.h>
#include <stdlib.h>

float uniformHost(float lb, float ub){
  float u = ((float) rand()) / ((float) RAND_MAX);
  return (ub - lb) * u + lb;
}