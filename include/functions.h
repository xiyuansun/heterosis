#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include<Chain.h>
#include<Config.h>

Config *config(int, char**);
void getopts(Config*, int, char**);
void printConfig(Config*);
void freeConfig(Config*);

float uniformHost(float, float);
float normalHost(float, float);
float gammaHost(float, float, float);
float betaHost(float, float);

#endif /* FUNCTIONS_H */