#include <functions.h>
#include <math.h>
#include <stdlib.h>

float normalHost(float m, float s){

  float u1 = uniformHost(0, 1);
  float u2 = uniformHost(0, 1);
  
  return sqrt(-2 * log(u1)) * sin(2 * M_PI * u2) * s + m;
}