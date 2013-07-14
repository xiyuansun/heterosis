#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <float.h>

float lfact(int n){
  return (n == 1 || n == 0) ? 1 : (lfact(n - 1) + (n ? log((float) n) : log(FLT_MIN)));
}

num_t logLik(count_t *y, int *group, int N, int G, 
              num_t *c, num_t *phi, num_t *alp, num_t *del, num_t *eps){
  int n, g;
  num_t ret = 0, mu, llam;
  
  for(n = 0; n < N; ++n){
    for(g = 0; g < G; ++g){
      
	  if(group[n] == 1){
		mu = phi[g] - alp[g];
	  } else if(group[n] == 2){
		mu = phi[g] + del[g];
	  } else {
		mu = phi[g] + alp[g];
	  }
	
      llam = c[n] + eps[iG(n, g)] + mu;
      ret += y[iG(n, g)] * llam - exp(llam) - lfact(y[iG(n, g)]);
    }
  }

  return ret;
}