#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void printParms(Chain *host_a, Chain *dev_a, Config *cfg){

  int m, n, g;
  int N = cfg->N, G = cfg->G;
  char file[BUF];
  FILE *fp;
  Chain *a;
  
  if(!cfg->parmsFlag)
    return;
      
  sprintf(file, "../out/parms/chain%d.txt", cfg->chainNum);
  fp = fopen(file, "w");
      
  if(fp == NULL){
    printf("ERROR: unable to create file, %s\n", file);
    return;
  } 
    
  for(n = 0; n < cfg->N; ++n)
    fprintf(fp, "c%d ", n);
    
  for(g = 0; g < cfg->G; ++g)
    fprintf(fp, "phi%d ", g);
    
  for(g = 0; g < cfg->G; ++g)
    fprintf(fp, "alpha%d ", g);
    
  for(g = 0; g < cfg->G; ++g)
    fprintf(fp, "delta%d ", g);
      
  for(g = 0; g < cfg->G; ++g)
    fprintf(fp, "eta%d ", g);
      
  for(g = 0; g < cfg->G; ++g)
    for(n = 0; n < cfg->N; ++n)
      fprintf(fp, "eps_lib%d_gene%d ", n, g);

  fprintf(fp, "\n");
  a = chainDeviceToHost(host_a, dev_a, cfg);
    
  for(m = 0; m <= cfg->M; ++m){
    for(n = 0; n < cfg->N; ++n)
      fprintf(fp, NUM_TF, a->c[iN(m, n)]); fprintf(fp, " ");
      
    for(g = 0; g < cfg->G; ++g)
      fprintf(fp, NUM_TF, phi[iG(m, g)]); fprintf(fp, " ");

    for(g = 0; g < cfg->G; ++g)
      fprintf(fp, NUM_TF, alp[iG(m, g)]); fprintf(fp, " ");
      
    for(g = 0; g < cfg->G; ++g)
      fprintf(fp, NUM_TF, del[iG(m, g)]); fprintf(fp, " ");

    for(g = 0; g < cfg->G; ++g)
      fprintf(fp, NUM_TF, eta[iG(m, g)]); fprintf(fp, " ");
      
    for(n = 0; n < cfg->N; ++n)
      for(g = 0; g < cfg->G; ++g)
        fprintf(fp, NUM_TF, eps[iNG(m, n, g)]); fprintf(fp, " ");
      
    fprintf(fp, "\n");
  } 
    
  fclose(fp);
}