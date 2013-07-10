#ifndef DEVICEFUNCTIONS_H
#define DEVICEFUNCTIONS_H

inline __device__ num_t mu(Chain *a, int n, num_t phi, num_t alp, num_t del){
  if(a->grp[n] == 1){
    return phi - alp;
  } else if(a->grp[n] == 2){
    return phi + del;
  } else {
    return phi + alp;
  }
}

inline __device__ num_t runiformDevice(Chain *a, int g, num_t lb, num_t ub){
  num_t u = curand_uniform(&(a->states[g]));
  return (ub - lb) * u + lb;
}

inline __device__ num_t rnormal(Chain *a, int g, num_t m, num_t s){

  num_t u1 = runiformDevice(a, g, 0, 1);
  num_t u2 = runiformDevice(a, g, 0, 1);
  
  return sqrt(-2 * log(u1)) * sin(2 * M_PI * u2) * s + m;
}

inline __device__ num_t rgammaDevice(Chain *a, int g, num_t shape, num_t rate, num_t lb){
   
  num_t A, c, d, r, u, v, w, x, z, eps, eps0, 
        haznaz, lam, ret, tmp1, tmp2;

  if(shape <= 0){
    return(-1);
  }
  
  if(rate <= 0){
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
      x = - (1/eps) * log(runiformDevice(a, g, 0, 1)) + lb * rate;
      u = runiformDevice(a, g, 0, 1);

      if(u < A * pow(x, shape - 1) * exp((eps - 1) * x))
        return(x / rate);
    }

  } else if(shape >= 1){ /* Marsaglia and Tsang (2000) */

    d = shape - 1/3;
    c = 1 / sqrt(9 * d);

    while(1){
      v = -1;
      while(v <= 0){
        x = rnormalDevice(a, g, 0, 1);
        v = pow(1 + c*x, 3);
      }

      ret = d * v / rate;

      if(ret > lb){
        u = runiformDevice(a, g, 0, 1);

        if(u < 1 - 0.0331 * pow(x, 4))
          return(ret);

        if(log(u) < 0.5 * pow(x, 2) + d * (1 - v + log(v)))
          return(ret);
      }
    }
  } else if (0.135 <= shape && shape < 1){ /* Kundu and Gupta (2006) */

    while(1){      

      u = runiformDevice(a, g, 0, 1);
      x = -2 * log(1 - pow(u, 1 / shape));
      ret = x / rate;

      if(ret > lb){
        v = runiformDevice(a, g, 0, 1);

        tmp1 = exp(-x/2);
        tmp2 = pow(x, shape - 1)* tmp1 * pow(((float) 2), ((float) (1 - shape))) * pow(((float) (1 - tmp1)), ((float) (1 - shape)));

        if(v < tmp2)
          return(ret);
      }
    }
  } else{ /* Martin and Liu (2013) */
   
    while(1){  
      lam = 1/shape - 1;
      w = shape / (exp(1 - shape));
      r = 1 / (1 + w); 
      u = runiformDevice(a, g, 0, 1);

      if(u <= r){
        z = - log(u / r);
      } else {
        z = log(runiformDevice(a, g, 0, 1)) / lam;
      }
      
      ret = exp(-z / shape) / rate;

      if(ret > lb){
        if(z >= 0){
          haznaz = exp(-exp(-z / shape));
        } else{
          haznaz = 1 / (w * lam) * exp((lam - 1) * z - exp(z / shape));
        }

        if(haznaz > runiformDevice(a, g, 0, 1))
          return(ret);
      }
    }
  }
}

inline __device__ num_t rbetaDevice(Chain *chain, int g, num_t a, num_t b){

  num_t x = rgammaDevice(chain, g, a, 1, 0);
  num_t y = rgammaDevice(chain, g, b, 1, 0);
  
  return x / (x + y);
}

#endif /* DEVICEFUNCTIONS_H */