#include <functions.h>
#include <constants.h>
#include <math.h>
#include <stdlib.h>

num_t rbeta(num_t a, num_t b){

  num_t x = rgamma(a, 1, 0);
  num_t y = rgamma(b, 1, 0);
  
  return x / (x + y);
}
