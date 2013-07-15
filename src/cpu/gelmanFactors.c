#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void gelmanFactors(Config *cfg){
  char file[BUF];
  FILE *fp;

  if(cfg->gelman){
  
    if(cfg->verbose)
      printf("Calculating Gelman potential scale reduction factors.\n");

    system("R CMD BATCH ../R/gelman-factors.r");
  }  
}