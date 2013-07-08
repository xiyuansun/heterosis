#include <constants.h>
#include <stdio.h>
#include <stdlib.h>

void pi1(int *v, int l, const char* m){
  int i;
  printf("%s", m);
  for(i = 0; i < l; ++i)
    printf("%d ", v[i]);
  printf("\n\n"); 
}

void pf1(num_t *v, int l, const char* m){
  int i;
  printf("%s", m);
  for(i = 0; i < l; ++i){
    printf(NUM_TF, v[i]);
    printf(" ");
  }
  printf("\n\n"); 
}

void pi2(int **v, int l1, int l2, const char* m){
  int i, j;
  printf("%s", m);
  for(i = 0; i < l1; ++i){
    for(j = 0; j < l2; ++j)
      printf("%d ", v[i][j]);
    printf("\n");
  }
  printf("\n"); 
}

void pf2(num_t **v, int l1, int l2, const char* m){
  int i, j;
  printf("%s", m);
  for(i = 0; i < l1; ++i){
    for(j = 0; j < l2; ++j){
      printf(NUM_TF, v[i][j]);
      printf(" ");
    }
    printf("\n");
  }
  printf("\n"); 
}

void pi3(int ***v, int l1, int l2, int l3, const char* m){
  int i, j, k;
  printf("%s", m);
  for(i = 0; i < l1; ++i){
    for(j = 0; j < l2; ++j){
      for(k = 0; k < l3; ++k)
        printf("%d ", v[i][j][k]);
      printf("\n");
    }
    printf("\n"); 
  }
}

void pf3(num_t ***v, int l1, int l2, int l3, const char* m){
  int i, j, k;
  printf("%s", m);
  for(i = 0; i < l1; ++i){
    for(j = 0; j < l2; ++j){
      for(k = 0; k < l3; ++k){
        printf(NUM_TF, v[i][j][k]);
        printf(" "); 
      }
      printf("\n");
    }
    printf("\n"); 
  }
}

