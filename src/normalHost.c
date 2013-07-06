#include <functions.h>
#include <math.h>
#include <numericTypes.h>
#include <stdlib.h>

num_t normalHost(num_t m, num_t s){

  num_t u1 = uniformHost(0, 1);
  num_t u2 = uniformHost(0, 1);
  
  return sqrt(-2 * log(u1)) * sin(2 * M_PI * u2) * s + m;
}