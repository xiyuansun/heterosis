#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <float.h>

#define BUF 256
#define MAXROW 16384
#define NUM_TF "%0.3f"
#define NUM_TMIN FLT_MIN

#define MAXTHREADS 512

#define G_GRID (cfg->G < MAXTHREADS ? cfg->G : MAXTHREADS)
#define G_BLOCK (((float) cfg->G) / NTHREADS)

#define N_THREAD (cfg->N < MAXTHREADS ? cfg->N : MAXTHREADS)
#define N_BLOCK (((float) cfg->N) / NTHREADS)

#define GN_GRID dim3(G_GRID, N_GRID)
#define GN_BLOCK dim3(G_BLOCK, N_BLOCK)

#define IDX ((blockDim.x * blockIdx.x) + threadIdx.x)
#define IDY ((blockDim.x * blockIdx.x) + threadIdx.x)

typedef int count_t;
typedef float num_t;

#endif /* CONSTANTS_H */