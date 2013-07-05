#include <functions.h>
#include <math.h>
#include <stdlib.h>

float betaHost(float a, float b){

  float x = gammaHost(a, 1, 0);
  float y = gammaHost(b, 1, 0);
  
  return x / (x + y);
}
