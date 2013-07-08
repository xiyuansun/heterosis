#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void runChain(Chain *a, Config *cfg){
  int m;
  
  for(m = 0; m < a->M; ++m){
    sampleC(a); printf("c\n");
    sampleTau(a, cfg); printf("sig-tau\n");
    samplePiAlp(a, cfg); printf("pialp\n");
    samplePiDel(a, cfg);printf("pidel\n");
    sampleD(a, cfg); printf("d\n");
    sampleThePhi(a, cfg); printf("thephi\n");
    sampleTheAlp(a, cfg); printf("thealp\n");
    sampleTheDel(a, cfg);printf("thedel\n");
    sampleSigC(a); printf("sigc\n");
    sampleSigPhi(a, cfg); printf("sigphi\n");
    sampleSigAlp(a, cfg);printf("sigalp\n");
    sampleSigDel(a, cfg);printf("sigdel\n");
    sampleEta(a); printf("eta\n");
    sampleEps(a); printf("eps\n");
    samplePhiAlpDel(a, cfg); printf("phialpdel\n");
  }
}