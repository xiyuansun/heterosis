#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <numericTypes.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv){
  int g, n, *s, **y;
  int *grp;
  
  Config *cfg = config(argc, argv);
  
  y = readData(cfg);
  
  grp = readGrp(cfg);
  

  s = malloc(5 * sizeof(int));
 
  
  for(n = 0; n < 10; ++n){
    s = mySampleIntHost(5, 5);
    
  /*  for(g = 0; g < 5; ++g)
      printf("%d ", s[g]);
    printf("\n"); */
  }  
  
  
  freeConfig(cfg);
  free(y);
  return 0;
}