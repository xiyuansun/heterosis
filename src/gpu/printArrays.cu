#include <constants.h>
#include <curand_kernel.h>
#include <stdio.h>
#include <stdlib.h>

__host__ void pi1(int *v, int l, const char* m){
  int i;
  printf("%s", m);
  for(i = 0; i < l; ++i)
    printf("%d ", v[i]);
  printf("\n\n"); 
}

__host__ void pf1(num_t *v, int l, const char* m){
  int i;
  printf("%s", m);
  for(i = 0; i < l; ++i){
    printf(NUM_TF, v[i]);
    printf(" ");
  }
  printf("\n\n"); 
}

__host__ void pi2(int *v, int l1, int l2, const char* m){
  int i, j;
  printf("%s", m);
  for(i = 0; i < l1; ++i){
    for(j = 0; j < l2; ++j)
      printf("%d ", v[i*l2 + j]);
    printf("\n");
  }
  printf("\n"); 
}

__host__ void pf2(num_t *v, int l1, int l2, const char* m){
  int i, j;
  printf("%s", m);
  for(i = 0; i < l1; ++i){
    for(j = 0; j < l2; ++j){
      printf(NUM_TF, v[i*l2 + j]);
      printf(" ");
    }
    printf("\n");
  }
  printf("\n"); 
}

__host__ void pstate(curandState *v, int l1, int l2, const char* m){
  int i, j;
  printf("%s", m);
  for(i = 0; i < l1; ++i){
    for(j = 0; j < l2; ++j)
      printf("%d ", (v[i*l2 + j]).d);
    printf("\n");
  }
  printf("\n"); 
}

__host__ void pi3(int *v, int l1, int l2, int l3, const char* m){
  int i, j, k;
  printf("%s", m);
  for(i = 0; i < l1; ++i){
    for(j = 0; j < l2; ++j){
      for(k = 0; k < l3; ++k)
        printf("%d ", v[i*l2*l3 + j*l3 + k]);
      printf("\n");
    }
    printf("\n"); 
  }
}

__host__ void pf3(num_t *v, int l1, int l2, int l3, const char* m){
  int i, j, k;
  printf("%s", m);
  for(i = 0; i < l1; ++i){
    for(j = 0; j < l2; ++j){
      for(k = 0; k < l3; ++k){
        printf(NUM_TF, v[i*l2*l3 + j*l3 + k]);
        printf(" "); 
      }
      printf("\n");
    }
    printf("\n"); 
  }
}

