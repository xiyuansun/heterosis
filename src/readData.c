#include <Config.h>
#include <constants.h>
#include <numericTypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

count_t **readData(Config *cfg){
  int g = 0, n = 0, **y;
  char *buf, row[MAXROW];
  FILE *fp = fopen(cfg->dataFile, "r");
  
  if(fp == NULL){
    printf("ERROR: data file \"%s\" not found.\n", cfg->dataFile);
    fclose(fp);
    return NULL;
  }

  if(!cfg->G || !cfg->N){
    while(fgets(row, MAXROW, fp) != NULL){ 
    
      if(!cfg->N && !n){
        buf = strtok(row, " \n");
      
        while(buf != NULL){    
          ++n;
          buf = strtok(NULL, " \n");
        }
        
        cfg->N = n;
      }
      
      if(cfg->G)
        break;
      
      ++g;
    }
    
    if(!cfg->G)
      cfg->G = g;
  }
  
  y = malloc(cfg->N * sizeof(count_t*));
  rewind(fp);
  
  for(n = 0; n < cfg->N; ++n)
    y[n] = malloc(cfg->G * sizeof(count_t)); 
  
  for(g = 0; g < cfg->G; ++g)
    for(n = 0; n < cfg->N; ++n)
      fscanf(fp, "%d", y[n] + g);
  
  fclose(fp);
  return y;
}