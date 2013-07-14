#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

__host__ void allocConfig(Config *cfg){
  cfg->yMeanG = (num_t*) malloc(cfg->N * sizeof(num_t));

  CUDA_CALL(cudaMalloc((void**) &(cfg->devY), cfg->N * cfg->G * sizeof(count_t)));
  CUDA_CALL(cudaMalloc((void**) &(cfg->devGrp), cfg->N * sizeof(count_t)));
  CUDA_CALL(cudaMalloc((void**) &(cfg->devYMeanG), cfg->N * sizeof(count_t)));
}