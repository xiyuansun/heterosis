#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

void updateDICprep(Chain *a){
  
  int n, g, N = a->N, G = a->G, m = a->m - a->burnin;
  num_t s1, s2, ll = logLik(a->y, a->grp, a->N, a->G, a->c, a->phi, a->alp, a->del, a->eps);
  
  if(m < 0){
    return;
  } else if(!m){
   
    for(n = 0; n < N; ++n)
      a->meanC[n] = a->c[n];
      
    for(g = 0; g < G; ++g){
      a->meanPhi[g] = a->phi[g];
      a->meanAlp[g] = a->alp[g];
      a->meanDel[g] = a->del[g];
      
      for(n = 0; n < N; ++n)
        a->meanEps[iG(n, g)] = a->eps[iG(n, g)];
    }
    
    a->meanLogLik = ll;
    
  } else {
    s1 = ((m - 1)/((num_t) m));
    s2 = (1/((num_t) m));
   
    for(n = 0; n < N; ++n)
      a->meanC[n] = s1 * a->meanC[n] + s2 * a->c[n];
      
    for(g = 0; g < G; ++g){
      a->meanPhi[g] = s1 * a->meanPhi[g] + s2 * a->phi[g];
      a->meanAlp[g] = s1 * a->meanAlp[g] + s2 * a->alp[g];
      a->meanDel[g] = s1 * a->meanDel[g] + s2 * a->del[g];
      
      for(n = 0; n < N; ++n)
        a->meanEps[iG(n, g)] = s1 * a->meanEps[iG(n, g)] + s2 * a->eps[iG(n, g)];
    }
    
    a->meanLogLik = s1 * (a->meanLogLik) + s2 * ll;
  }
}