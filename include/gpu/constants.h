#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <float.h>

#define BUF 256
#define MAXROW 16384
#define NUM_TF "%0.3f"
#define NUM_TMIN FLT_MIN

#define NTHREADS (cfg->G < cfg->nthreads ? cfg->G : cfg->nthreads)
#define NBLOCKS ceil(cfg->G / NTHREADS)
#define GENE ((blockDim.x * blockIdx.x) + threadIdx.x)

typedef int count_t;
typedef float num_t;

#endif /* CONSTANTS_H */