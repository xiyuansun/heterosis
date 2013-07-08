#include <functions.h>
#include <constants.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

num_t gammaHost(num_t shape, num_t rate, num_t lb){
   
  num_t A, c, d, r, u, v, w, x, z, eps, eps0, 
        haznaz, lam, ret, tmp1, tmp2;

  if(shape <= 0){
    printf("Error: bad shape: ");
    printf(NUM_TF, shape);
    printf("\n");
    return(-1);
  }
  
  if(rate <= 0){
    printf("Error: bad rate: ");
    printf(NUM_TF, rate);
    printf("\n");
    return(-1);
  }

  if((shape >= 1) && (shape - 1 < lb * rate)){ /* Chung (1998) */

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
        u = uniformHost(0, 1);

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
        v = uniformHost(0, 1);

        tmp1 = exp(-x/2);
        tmp2 = pow(x, shape - 1)* tmp1 * pow(2, 1 - shape) * pow(1 - tmp1, 1 - shape);

        if(v < tmp2)
          return(ret);
      }
    }
  } else{ /* Martin and Liu (2013) */
   
    while(1){  
      lam = 1/shape - 1;
      w = shape / (exp(1 - shape));
      r = 1 / (1 + w); 
      u = uniformHost(0, 1);

      if(u <= r){
        z = - log(u / r);
      } else {
        z = log(uniformHost(0, 1)) / lam;
      }
      
      ret = exp(-z / shape) / rate;

      if(ret > lb){
        if(z >= 0){
          haznaz = exp(-exp(-z / shape));
        } else{
          haznaz = 1 / (w * lam) * exp((lam - 1) * z - exp(z / shape));
        }

        if(haznaz > uniformHost(0, 1))
          return(ret);
      }
    }
  }
}