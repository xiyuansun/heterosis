#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

void dic(Chain *a){
  a->logLikMean = logLik(a->y, a->grp, a->N, a->G, a->meanC, a->meanPhi, 
                         a->meanAlp, a->meanDel, a->meanEps);
  a->dic = 4 * a->meanLogLik - 2 * a->logLikMean;
}