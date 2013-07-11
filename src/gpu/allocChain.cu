#include <Chain.h>
#include <Config.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

__host__ void allocChainHost(Chain **a, Config *cfg){

  *a = (Chain*) malloc(sizeof(Chain));
  fprintf(cfg->log, "  Allocating chain.\n"); 

  /* data */  
  
  (*a)->y = (count_t*) malloc(cfg->N * cfg->G * sizeof(count_t));
  (*a)->yMeanG = (num_t*) malloc(cfg->N * sizeof(num_t));
  (*a)->grp = (int*) malloc(cfg->N * sizeof(int));

  /* curand states */
  
  (*a)->states = (curandState*) malloc(cfg->G * sizeof(curandState));

  /* parameters */

  (*a)->c      = (num_t*) malloc((cfg->M + 1) * cfg->N * sizeof(num_t));
  (*a)->sigC   = (num_t*) malloc((cfg->M + 1) * sizeof(num_t));
  (*a)->eps    = (num_t*) malloc((cfg->M + 1) * cfg->N * cfg->G * sizeof(num_t));
  (*a)->eta    = (num_t*) malloc((cfg->M + 1) * cfg->G * sizeof(num_t));
  (*a)->d      = (num_t*) malloc((cfg->M + 1) * sizeof(num_t));  
  (*a)->tau    = (num_t*) malloc((cfg->M + 1) * sizeof(num_t));
  (*a)->phi    = (num_t*) malloc((cfg->M + 1) * cfg->G * sizeof(num_t));
  (*a)->thePhi = (num_t*) malloc((cfg->M + 1) * sizeof(num_t));
  (*a)->sigPhi = (num_t*) malloc((cfg->M + 1) * sizeof(num_t));
  (*a)->alp    = (num_t*) malloc((cfg->M + 1) * cfg->G * sizeof(num_t));
  (*a)->theAlp = (num_t*) malloc((cfg->M + 1) * sizeof(num_t));
  (*a)->sigAlp = (num_t*) malloc((cfg->M + 1) * sizeof(num_t));
  (*a)->piAlp  = (num_t*) malloc((cfg->M + 1) * sizeof(num_t));
  (*a)->del    = (num_t*) malloc((cfg->M + 1) * cfg->G * sizeof(num_t));
  (*a)->theDel = (num_t*) malloc((cfg->M + 1) * sizeof(num_t));
  (*a)->sigDel = (num_t*) malloc((cfg->M + 1) * sizeof(num_t));
  (*a)->piDel  = (num_t*) malloc((cfg->M + 1) * sizeof(num_t));
  
  /* temporary and return values */
  
  (*a)->tmp1 = (num_t*) malloc(cfg->G * sizeof(num_t));
  (*a)->tmp2 = (num_t*) malloc(cfg->G * sizeof(num_t));

  (*a)->Old   = (num_t*) malloc(cfg->N * sizeof(num_t));
  (*a)->New   = (num_t*) malloc(cfg->N * sizeof(num_t));
  (*a)->lOld  = (num_t*) malloc(cfg->N * sizeof(num_t));
  (*a)->lNew  = (num_t*) malloc(cfg->N * sizeof(num_t));

  /* tuning parameters for Metropolis steps */
  
  (*a)->tuneC = (num_t*) malloc(cfg->N * sizeof(num_t));
  (*a)->tunePhi = (num_t*) malloc(cfg->G * sizeof(num_t));
  (*a)->tuneEps = (num_t*) malloc(cfg->N * cfg->G * sizeof(num_t));

  /* number of acceptances for Metropolis steps */

  (*a)->accC   = (int*) malloc(cfg->N * sizeof(int));
  (*a)->accPhi = (int*) malloc(cfg->G * sizeof(int));
  (*a)->accAlp = (int*) malloc(cfg->G * sizeof(int));
  (*a)->accDel = (int*) malloc(cfg->G * sizeof(int));
  (*a)->accEps = (int*) malloc(cfg->N * cfg->G * sizeof(int));

}

__host__ void allocChainDevice(Chain **host_a, Chain **dev_a, Config *cfg){

  float myTime;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  *host_a = (Chain*) malloc(sizeof(Chain));
  
  /* data */  
    
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->y), cfg->N * cfg->G * sizeof(count_t))); 
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->yMeanG), cfg->N * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->grp), cfg->N * sizeof(int)));

  /* curand states */
  
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->states), cfg->N * cfg->G * sizeof(curandState)));

  /* parameters */

  CUDA_CALL(cudaMalloc((void**) &((*host_a)->c), (cfg->M + 1) * cfg->N * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->sigC), (cfg->M + 1) * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->eps), (cfg->M + 1) * cfg->N * cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->eta), (cfg->M + 1) * cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->d), (cfg->M + 1) * sizeof(num_t)));  
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->tau), (cfg->M + 1) * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->phi), (cfg->M + 1) * cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->thePhi), (cfg->M + 1) * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->sigPhi), (cfg->M + 1) * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->alp), (cfg->M + 1) * cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->theAlp), (cfg->M + 1) * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->sigAlp), (cfg->M + 1) * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->piAlp), (cfg->M + 1) * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->del), (cfg->M + 1) * cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->theDel), (cfg->M + 1) * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->sigDel), (cfg->M + 1) * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->piDel), (cfg->M + 1) * sizeof(num_t)));
  
  /* temporary and return values */
  
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->tmp1), cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->tmp2), cfg->G * sizeof(num_t)));

  CUDA_CALL(cudaMalloc((void**) &((*host_a)->Old), cfg->N * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->New), cfg->N * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->lOld), cfg->N * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->lNew), cfg->N * sizeof(num_t)));

  /* tuning parameters for Metropolis steps */
  
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->tuneC), cfg->N * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->tunePhi), cfg->G * sizeof(num_t)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->tuneEps), cfg->N * cfg->G * sizeof(num_t)));

  /* number of acceptances for Metropolis steps */

  CUDA_CALL(cudaMalloc((void**) &((*host_a)->accC), cfg->N * sizeof(int)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->accPhi), cfg->G * sizeof(int)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->accAlp), cfg->G * sizeof(int)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->accDel), cfg->G * sizeof(int)));
  CUDA_CALL(cudaMalloc((void**) &((*host_a)->accEps), cfg->N * cfg->G * sizeof(int)));
    
  CUDA_CALL(cudaMalloc((void**) dev_a, sizeof(Chain)));  
  CUDA_CALL(cudaMemcpy(*dev_a, *host_a, sizeof(Chain), cudaMemcpyHostToDevice));
  
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  
  fprintf(cfg->time, "%0.3f ", myTime); /* elapsed time in minutes */
}