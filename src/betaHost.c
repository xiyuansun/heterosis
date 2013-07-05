#include <functions.h>
#include <math.h>
#include <stdlib.h>

float betaHost(float a, float b){

  float x = sampleGammaHost(a, 1, 0);
  float y = sampleGammaHost(b, 1, 0);
  
  return x / (x + y);
}
