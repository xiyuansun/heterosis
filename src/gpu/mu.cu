#include <Chain.h>
#include <constants.h>
#include <stdlib.h>

inline __device__ num_t mu(Chain *a, int n, num_t phi, num_t alp, num_t del){
  if(a->grp[n] == 1){
    return phi - alp;
  } else if(a->grp[n] == 2){
    return phi + del;
  } else {
    return phi + alp;
  }
}
