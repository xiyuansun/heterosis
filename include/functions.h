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

num_t lEps(Chain*, int, int, num_t);
void sampleEps_kernel1(Chain*);
void sampleEps_kernel2(Chain*);
void sampleEps(Chain*);

void lD_kernel1(Chain*);
void lD_kernel2(Chain*);
void lD_kernel3(Chain*);
void lD_kernel4(Chain*, int);
void lD(Chain*, int);
void sampleD_kernel1(Chain*);
void sampleD_kernel2(Chain*);
void sampleD(Chain*, Config*);

num_t lPhi(Chain*, int, num_t);
void samplePhi_kernel1(Chain*);
void samplePhi_kernel2(Chain*);
void samplePhi(Chain*);

num_t alpProp(Chain*, int);
num_t lAlp(Chain*, int, num_t);
void sampleAlp_kernel1(Chain*);
void sampleAlp_kernel2(Chain*);
void sampleAlp(Chain*);

num_t delProp(Chain*, int);
num_t lDel(Chain*, int, num_t);
void sampleDel_kernel1(Chain*);
void sampleDel_kernel2(Chain*);
void sampleDel(Chain*);

num_t lPhiAlpDelJoint(Chain*, int, num_t, num_t, num_t);
void samplePhiAlpDel_kernel1Joint(Chain*);
void samplePhiAlpDel_kernel2Joint(Chain*);
void samplePhiAlpDelJoint(Chain*);

void samplePhiAlpDel(Chain*, Config*);

void sampleSigC(Chain*);

void sampleEta_kernel1(Chain*);
void sampleEta_kernel2(Chain*);
void sampleEta_kernel3(Chain*);
void sampleEta(Chain*);

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
void oneChain(int, char**);

void printProbs(Chain*, Config*);
void printRates(Chain*, Config*);
void printHyper(Chain*, Config*);
void printParms(Chain*, Config*);
void summarizeChain(Chain*, Config*);

#endif /* FUNCTIONS_H */