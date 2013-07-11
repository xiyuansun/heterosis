#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <float.h>

#define BUF 256
#define MAXROW 16384
#define NUM_TF "%0.3f"
#define NUM_TMIN FLT_MIN

#define MAXTHREADS 512
#define NTHREADS (cfg->G < MAXTHREADS ? cfg->G : MAXTHREADS)
#define NBLOCKS ceil(((float) cfg->G) / NTHREADS)
#define GENE ((blockDim.x * blockIdx.x) + threadIdx.x)

typedef int count_t;
typedef float num_t;

#endif /* CONSTANTS_H */