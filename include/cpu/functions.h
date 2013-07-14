#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "Chain.h"
#include "Config.h"
#include "constants.h"

#define iG(n, g) ((n) * G + (g))

void pi1(int*, int, const char*);
void pf1(num_t*, int, const char*);
void pi2(int*, int, int, const char*);
void pf2(num_t*, int, int, const char*);
void pi3(int*, int, int, int, const char*);
void pf3(num_t*, int, int, int, const char*);

Config *config(int, char**);
void getopts(Config*, int, char**);
void printConfig(Config*);

int *mySampleInt(int, int);
int *readGrp(Config*);
count_t *readData(Config*);

Chain *allocChain(Config*);
Chain *newChain(Config*);
void newChain_kernel1(Chain*);
void newChain_kernel2(Chain*);
void printChain(Chain*);
void freeChain(Chain*, Config*);

num_t mu(Chain*, int, num_t, num_t, num_t);
num_t runiform(num_t, num_t);
num_t rnormal(num_t, num_t);
num_t rgamma(num_t, num_t, num_t);
num_t rbeta(num_t, num_t);

void lC_kernel1(Chain*, int);
void lC_kernel2(Chain*, int);
void lC_kernel3(Chain*, int, int);
void lC(Chain*, int, int);
void sampleC_kernel1(Chain*);
void sampleC_kernel2(Chain*);
void sampleC(Chain*, Config*);

num_t lEps(Chain*, int, int, num_t);
void sampleEps_kernel(Chain*);
void sampleEps(Chain*, Config*);

void lD_kernel1(Chain*);
void lD_kernel2(Chain*);
void lD_kernel3(Chain*);
void lD_kernel4(Chain*, int);
void lD(Chain*, int);
void sampleD_kernel1(Chain*);
void sampleD_kernel2(Chain*);
void sampleD(Chain*, Config*);

num_t lPhi(Chain*, int, num_t);
void samplePhi_kernel(Chain*);
void samplePhi(Chain*, Config*);

num_t alpProp(Chain*, int);
num_t lAlp(Chain*, int, num_t);
void sampleAlp_kernel(Chain*);
void sampleAlp(Chain*, Config*);

num_t delProp(Chain*, int);
num_t lDel(Chain*, int, num_t);
void sampleDel_kernel(Chain*);
void sampleDel(Chain*, Config*);

num_t lPhiAlpDelJoint(Chain*, int, num_t, num_t, num_t);
void samplePhiAlpDelJoint_kernel(Chain*);
void samplePhiAlpDelJoint(Chain*, Config*);

void samplePhiAlpDel(Chain*, Config*);

void sampleSigC(Chain*, Config*);

void sampleEta_kernel1(Chain*);
void sampleEta_kernel2(Chain*);
void sampleEta(Chain*, Config*);

void sampleTau_kernel1(Chain*);
void sampleTau_kernel2(Chain*);
void sampleTau_kernel3(Chain*);
void sampleTau(Chain*, Config*);

void sampleThePhi_kernel1(Chain*);
void sampleThePhi_kernel2(Chain*);
void sampleThePhi(Chain*, Config*);

void sampleTheAlp_kernel1(Chain*);
void sampleTheAlp_kernel2(Chain*);
void sampleTheAlp_kernel3(Chain*);
void sampleTheAlp_kernel4(Chain*);
void sampleTheAlp(Chain*, Config*);

void sampleTheDel_kernel1(Chain*);
void sampleTheDel_kernel2(Chain*);
void sampleTheDel_kernel3(Chain*);
void sampleTheDel_kernel4(Chain*);
void sampleTheDel(Chain*, Config*);

void sampleSigPhi_kernel1(Chain*);
void sampleSigPhi_kernel2(Chain*);
void sampleSigPhi_kernel3(Chain*);
void sampleSigPhi(Chain*, Config*);

void sampleSigAlp_kernel1(Chain*);
void sampleSigAlp_kernel2(Chain*);
void sampleSigAlp_kernel3(Chain*);
void sampleSigAlp_kernel4(Chain*);
void sampleSigAlp(Chain*, Config*);

void sampleSigDel_kernel1(Chain*);
void sampleSigDel_kernel2(Chain*);
void sampleSigDel_kernel3(Chain*);
void sampleSigDel_kernel4(Chain*);
void sampleSigDel(Chain*, Config*);

void samplePiAlp_kernel1(Chain*);
void samplePiAlp_kernel2(Chain*);
void samplePiAlp_kernel3(Chain*);
void samplePiAlp(Chain*, Config*);

void samplePiDel_kernel1(Chain*);
void samplePiDel_kernel2(Chain*);
void samplePiDel_kernel3(Chain*);
void samplePiDel(Chain*, Config*);

void runChain(Chain*, Config*);
void oneChain(Config*);
void chains(int, char**);

void printHeaders(Chain*, Config*);
void sumLogLik_kernel(Chain*);
void intermResults(Chain*, Config*);
void summarizeChain(Chain*, Config*);

void logLiks_kernel1(Chain*);
void logLiks_kernel2(Chain*);
void logLiks(Chain*, Config*);

#endif /* FUNCTIONS_H */
