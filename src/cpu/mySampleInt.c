#include <stdio.h>
#include <stdlib.h>

int *mySampleInt(int N, int n){

  int i, j, repeats, *ret = malloc(n * sizeof(int)), swap, tmp;

  /* sample without replacement */

  for(i = 0; i < n; ++i)
    ret[i] = -1;

  i = 0;
  while(i < n){
    ret[i] = rand() % N;

    repeats = 0;
    j = 1;

    for(j = 0; j < i; ++j)
      if(ret[j] == ret[i])
        ++repeats;

    if(!repeats)
      ++i;
  }

  /* bubble sort */

  if(n < 2)
    return(ret);

  swap = 1;
  while(swap){
    swap = 0;
    
    for(i = 1; i < n; ++i){
      if(ret[i - 1] > ret[i]){
        tmp = ret[i]; 
        ret[i] = ret[i - 1];
        ret[i - 1] = tmp;
        
        swap = 1;
      }
    }
  }

  return ret;
}