#ifndef CONSTANTS_H
#define CONSTANTS_H

#define USE_R 0 /* set to 1 to use R script, gelman-factors.r, 
                 * to compute Gelman factors.
                 */

#define BUF 256
#define MAXROW 16384
#define NUM_TMIN (-1.0/0.0)
#define M_PI 3.1415926535897932384626433832795028841971693993751058
#define MAX_NG ((N > G) ? N : G)
#define SECONDS 0.001 /* number of seconds in one unit of output time */

typedef int count_t;

typedef double num_t;
#define NUM_TF "%f"


#endif /* CONSTANTS_H */