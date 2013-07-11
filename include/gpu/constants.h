#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <float.h>

#define BUF 256
#define MAXROW 16384
#define NUM_TF "%0.3f"
#define NUM_TMIN FLT_MIN

#define MAXTHREADS 512

#define G_BLOCK (cfg->G < MAXTHREADS ? cfg->G : MAXTHREADS)
#define G_GRID ((int) ceil(((float) cfg->G) / G_BLOCK))

#define N_BLOCK (cfg->N < MAXTHREADS ? cfg->N : MAXTHREADS)
#define N_GRID ((int) ceil(((float) cfg->N) / N_BLOCK))

#define GN_GRID dim3(G_GRID, N_GRID, 1)
#define GN_BLOCK dim3(G_BLOCK, N_BLOCK, 1)

#define IDX ((blockDim.x * blockIdx.x) + threadIdx.x)
#define IDY ((blockDim.x * blockIdx.x) + threadIdx.x)
#define ID ((gridDim.x * blockDim.x * IDY) + IDX)

typedef int count_t;
typedef float num_t;

#endif /* CONSTANTS_H */ 