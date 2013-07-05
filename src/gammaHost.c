#include <functions.h>
#include <math.h>
#include <stdlib.h>

float gammaHost(float shape, float rate, float lb){
   
  float A, c, d, u, v, w, x, z, eps, eps0, lam;

  if(shape <= 0){
    printf("Error: bad shape: %0.3f", shape);
    return(0/0);
  }
  
  if(rate <= 0){
    printf("Error: bad rate: %0.3f", rate);
    return(0/0);
  }

  if(shape - 1 < lb * rate){ /* Chung (1998) */

    c = lb * rate;
    eps0 = (c - shape + sqrt(pow(c - shape, 2) + 4 * c))/(2 * c);

    if(eps0 > 1){
      eps = 0.75;
    } else {
      eps = eps0;
    }

    A = pow((1 - eps)/(shape - 1), shape - 1) * exp(shape - 1);

    while(1){
      x = - (1/eps) * log(uniformHost(0, 1)) + lb * rate;
      u = uniformHost(0, 1);

      if(u < A * pow(x, shape - 1) * exp((eps - 1) * x))
        return(x / rate);
    }

  } else if(shape >= 1){ /* Marsaglia and Tsang (2000) */

    d = shape - 1/3;
    c = 1 / sqrt(9 * d);

    while(1){
      v = -1;
      while(v <= 0){
        x = normalHost(0, 1);
        v = pow(1 + c*x, 3);
      }

      ret = d * v / rate;

      if(ret > lb){
        u = uniformHost(1);

        if(u < 1 - 0.0331 * pow(x, 4))
          return(ret);

        if(log(u) < 0.5 * pow(x, 2) + d * (1 - v + log(v)))
          return(ret);
      }
    }
  } else if (0.135 <= shape && shape < 1){ /* Kundu and Gupta (2006) */

    while(1){      

      u = uniformHost(0, 1);
      x = -2 * log(1 - pow(u, 1 / shape));
      ret = x / rate;

      if(ret > lb){
        v = runiformHost(1);

        tmp1 = exp(-x/2);
        tmp2 = pow(x, shape - 1)* tmp1 * pow(2, 1 - shape) * pow(1 - tmp1, 1 - shape);

        if(v < tmp2)
          return(ret);
      }
    }
  } else{ /* Martin and Liu (2013) */
   
    while(1){ # 
      lam = 1/shape - 1;
      w = shape / (exp(1 - shape));
      r = 1 / (1 + w); 
      u = sampleUniform(0, 1);

      if(u <= r){
        z = - log(u / r);
      } else {
        z = log(sampleUniform(0, 1)) / lam;
      }
      
      ret = exp(-z / shape) / rate;

      if(ret > lb){
        if(z >= 0){
          haznaz = exp(-exp(-z / shape));
        } else{
          haznaz = 1 / (w * lam) * exp((lam - 1) * z - exp(z / shape));
        }

        if(haznaz > sampleUniform(0, 1))
          return(ret);
      }
    }
  }
}