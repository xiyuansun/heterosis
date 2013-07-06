#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <Chain.h>
#include <Config.h>
#include <numericTypes.h>

void pi1(int*, int, const char*);
void pf1(float*, int, const char*);
void pi2(int**, int, int, const char*);
void pf2(float**, int, int, const char*);
void pi3(int***, int, int, int, const char*);
void pf3(float***, int, int, int, const char*);

Config *config(int, char**);
void getopts(Config*, int, char**);
void printConfig(Config*);
void freeConfig(Config*);

int *mySampleIntHost(int, int);
int *readGrp(Config*);
count_t **readData(Config*);

Chain *allocChainHost(Config*);
void newChainHost_kernel1(Chain*);
void newChainHost_kernel2(Chain*);
Chain *newChainHost(Config*);
void printChain(Chain*);
void freeChainHost(Chain*, Config*);

num_t uniformHost(num_t, num_t);
num_t normalHost(num_t, num_t);
num_t gammaHost(num_t, num_t, num_t);
num_t betaHost(num_t, num_t);

#endif /* FUNCTIONS_H */