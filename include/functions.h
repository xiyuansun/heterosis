#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <Chain.h>
#include <Config.h>
#include <constants.h>

void pi1(int*, int, const char*);
void pf1(num_t*, int, const char*);
void pi2(int**, int, int, const char*);
void pf2(num_t**, int, int, const char*);
void pi3(int***, int, int, int, const char*);
void pf3(num_t***, int, int, int, const char*);

Config *config(int, char**);
void getopts(Config*, int, char**);
void printConfig(Config*);
void freeConfig(Config*);

int *mySampleIntHost(int, int);
int *readGrp(Config*);
count_t **readData(Config*);

Chain *allocChainHost(Config*);
Chain *newChainHost(Config*);
void newChainHost_kernel1(Chain*);
void newChainHost_kernel2(Chain*);
void printChain(Chain*);
void freeChainHost(Chain*, Config*);

num_t mu(Chain*, int, num_t, num_t, num_t);
num_t uniformHost(num_t, num_t);
num_t normalHost(num_t, num_t);
num_t gammaHost(num_t, num_t, num_t);
num_t betaHost(num_t, num_t);

void lC_kernel1(Chain*, int);
void lC_kernel2(Chain*, int);
void lC_kernel3(Chain*, int, int);
void lC(Chain*, int, int);
void sampleC_kernel1(Chain*);
void sampleC_kernel2(Chain*);
void sampleC_kernel3(Chain*);
void sampleC(Chain*);

#endif /* FUNCTIONS_H */