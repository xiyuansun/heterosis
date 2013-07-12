#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

count_t *readData(Config *cfg){
  int g = 0, n = 0, G, *y;
  char *buf, row[MAXROW];

  FILE *fp = fopen(cfg->dataFile, "r");
  
  if(fp == NULL){
    printf("ERROR: data file \"%s\" not found.\n", cfg->dataFile);
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
  
  y = (count_t*) malloc(cfg->N * cfg->G * sizeof(count_t));
  rewind(fp);
  
  G = cfg->G;
  
  for(g = 0; g < cfg->G; ++g)
    for(n = 0; n < cfg->N; ++n)
      fscanf(fp, "%d", y + iG(n, g));
      
  fclose(fp);
  return y;
}