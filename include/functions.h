#ifndef FUNCTIONS_H

#include<Chain.h>
#include<Config.h>
#include<constants.h>

#define

Config *config(int, char**);
void getopts(Config*, int, char**);
void printConfig(Config*);
void freeConfig(Config*);

float uniformHOst(float, float);
float normalHost(float, float);
float gammaHost(float, float, float);
float betaHost(float, float);

#endif /* FUNCTIONS_H */