#include <Chain.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <thrust/reduce.h>

__global__ void lC_kernel1(Chain *a, int n){ /* kernel <<<G, 1>>> */
  int g = IDX, G = a->G;
  
  if(g < G)
    a->tmp1[g] = exp(a->eps[iG(n, g)] + mu(a, n, a->phi[g], 
                    a->alp[g], a->del[g]));
}

__global__ void lC_kernel2(Chain *a, int n, int newArg){ /* kernel <<<1, 1>>> */
  num_t arg, ret;

  if(newArg){
    arg = a->New[n];
  } else {
    arg = a->Old[n];
  }

  ret = arg * a->G * a->yMeanG[n] - exp(arg) * a->s1 - (arg*arg) / 
        (2 * a->sigC * a->sigC);

  if(newArg){
    a->lNew[n] = ret;
  } else {
    a->lOld[n] = ret;
  }
}

__host__ void lC(Chain *host_a, Chain *dev_a, Config *cfg, int n, int newArg){ /* host */
  lC_kernel1<<<G_GRID, G_BLOCK>>>(dev_a, n);
  
  thrust::device_ptr<num_t> tmp1(host_a->tmp1);  
  num_t s1 = thrust::reduce(tmp1, tmp1 + cfg->G);
  CUDA_CALL(cudaMemcpy(&(dev_a->s1), &s1, sizeof(num_t), cudaMemcpyHostToDevice));
  
  lC_kernel2<<<1, 1>>>(dev_a, n, newArg);
}

__global__ void sampleC_kernel1(Chain *a){ /* kernel <<<1, N>>> */
  int n = IDX;
  
  if(n < a->N){
    a->Old[n] = a->c[n];
    a->New[n] = rnormalDevice(a, n, a->Old[n], a->tuneC[n]);
    
  }
}

__global__ void sampleC_kernel2(Chain *a){ /* kernel <<<1, N>>> */
  int n = IDX;
  num_t dl, lp, lu;

  if(n < a->N){ 

    dl = a->lNew[n] - a->lOld[n];
    lp = 0 < dl ? 0 : dl;
    lu = log(runiformDevice(a, n, 0, 1));
      
    if(lu < lp){ /* accept */
      a->c[n] = a->New[n];
      a->tuneC[n] *= 1.1; /* Increase the proposal variance to avoid  
                                       gettiG stuck in a mode */
                                       
      if(a->m > a->burnin)                                 
        ++a->accC[n];
    } else { /* reject */
      a->tuneC[n] /= 1.1; /* If you're rejectiG too often, decrease the  
                                       proposal variance to sample closer to 
                                       the last accepted value. */
    }
    
    if(a->m > a->burnin) 
      a->sumC[n] += a->c[n]; 
  }
}

__host__ void sampleC(Chain *host_a, Chain *dev_a, Config *cfg){ /* host */
  int n;
  num_t myTime;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);  
  
  if(cfg->verbose)
    printf("c ");
   
  sampleC_kernel1<<<N_GRID, N_BLOCK>>>(dev_a);

  for(n = 0; n < cfg->N; ++n){ 
    lC(host_a, dev_a, cfg, n, 1);
    lC(host_a, dev_a, cfg, n, 0);
  }

  sampleC_kernel2<<<N_GRID, N_BLOCK>>>(dev_a);

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&myTime, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  
  cfg->timeC = myTime / MILLISECS;
}