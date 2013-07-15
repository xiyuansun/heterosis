#ifndef CONSTANTS_H
#define CONSTANTS_H

#define BUF 256
#define MAXROW 16384
#define NUM_TMIN (-1.0/0.0)
#define MAX_NG ((N > G) ? N : G)

#define MAXTHREADS 512
#define IDX ((blockDim.x * blockIdx.x) + threadIdx.x)

#define N_BLOCK (cfg->N < MAXTHREADS ? cfg->N : MAXTHREADS)
#define N_GRID ((cfg->N / N_BLOCK) + 1)

#define G_BLOCK (cfg->G < MAXTHREADS ? cfg->G : MAXTHREADS)
#define G_GRID ((cfg->G / G_BLOCK) + 1)

#define NG_BLOCK (MAX_NG < MAXTHREADS ? MAX_NG : MAXTHREADS)
#define NG_GRID ((MAX_NG / NG_BLOCK) + 1)

#define MILLISECS 1.0 /* number of milliseconds in a unit of time in output */

typedef int count_t;
typedef double num_t;
#define NUM_TF "%0.3f"

#endif /* CONSTANTS_H */