#include <functions.h>
#include <math.h>
#include <numericTypes.h>
#include <stdlib.h>

num_t betaHost(num_t a, num_t b){

  num_t x = gammaHost(a, 1, 0);
  num_t y = gammaHost(b, 1, 0);
  
  return x / (x + y);
}
