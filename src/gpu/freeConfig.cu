#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

__host__ void allocConfig(Config *cfg){
  free(cfg->y);
  free(cfg->grp);
  free(cfg->yMeanG);
  
  cudaFree(cfg->devY);
  cudaFree(cfg->devGrp);
  cudaFree(cfg->devYMeanG);

  free(cfg);  
}