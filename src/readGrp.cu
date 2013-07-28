#include <Config.h>
#include <stdio.h>
#include <stdlib.h>

int *readGrp(Config *cfg){

  int i, j, n, *grp, match = 0, nunique = 0, *unique;
  FILE *fp;

  if(cfg->N < 2){
    fprintf(stderr, "ERROR: bad experimental design.\n");
    exit(EXIT_FAILURE);
  }
  
  fp = fopen(cfg->groupFile, "r");
  
  if(fp == NULL){
    fprintf(stderr, "ERROR: group file \"%s\" not found.\n", cfg->groupFile);
    return NULL;
  }
  
  grp = (int*) malloc(cfg->N * sizeof(int));
  unique = (int*) calloc(cfg->N, sizeof(int));
  
  for(n = 0; n < cfg->N; ++n)
    fscanf(fp, "%d", grp + n);

  for(i = 0; i < cfg->N; ++i){ 
    match = 0;

    for(j = 0; j < i; ++j)
      if(grp[i] == unique[j])
        ++match;

    if(!match){
      ++nunique;
      unique[nunique - 1] = grp[i];
    } 
  }

  if (nunique == 2){

    cfg->heterosis = 0;
    
    for(n = 0; n < cfg->N; ++n){
      if(grp[n] == unique[1]){
        grp[n] = 1;
      } else{
        grp[n] = 3;
      }
    }
  } else if (nunique == 3){
    cfg->heterosis = 1;
    
  } else {  
    fprintf(stderr, "ERROR: bad experimental design.");
    fclose(fp);
    free(grp);
    free(unique);
    exit(EXIT_FAILURE);
  }
  
  fclose(fp);
  free(unique);
  return grp;
}
