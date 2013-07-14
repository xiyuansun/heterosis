#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdlib.h>

num_t rnormal(num_t m, num_t s){

  num_t u1 = runiform(0, 1);
  num_t u2 = runiform(0, 1);
  
  return sqrt(-2 * log(u1)) * sin(2 * M_PI * u2) * s + m;
}