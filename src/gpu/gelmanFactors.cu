#include <constants.h>
#include <functions.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

__host__ int getNumFiles(char *dir){
  FILE *fp;
  char buf[MAXROW], cmd[BUF];
  
  sprintf(cmd, "ls %s | wc -l", dir);
  fp = popen(cmd, "r");
  fgets(buf, MAXROW, fp);
  
  pclose(fp);
  return atoi(buf);
}

__host__ void extractOneParm(num_t **x, char *dir, FILE *outfp, char *parmName, 
                    int parmNum,  int nRows, int J, int burnin){
  FILE *fp;
  char cmd[BUF], buf[BUF];
  int i, j;
  
  for(j = 0; j < J; ++j){
    sprintf(cmd, "cut -d \' \' -f %d %schain%d.txt", parmNum, dir, j);
    /* sprintf(cmd, "../sh/merge-cols.sh %s %d", dir, j); */ /* <-- much slower */
    
    fp = popen(cmd, "r");
    
    for(i = 0; i < nRows; ++i){
      fgets(buf, BUF, fp);
      if(i > burnin)
        x[i - burnin - 1][j] = atof(buf); 
    }
    
    pclose(fp);
  }
}

__host__ num_t gelmanFactor(num_t **x, int I, int J){
  int i, j;
  num_t tmp, x__ = 0, *x_j, W = 0, B = 0;  
  
  x_j = (num_t*) calloc(J, sizeof(num_t));   
  
  for(i = 0; i < I; ++i){  
    for(j = 0; j < J; ++j){
      tmp = x[i][j];
      x__ += tmp;
      x_j[j] += tmp;
    }
  }
  
  x__ /= I * J; 
    
  for(j = 0; j < J; ++j)
    x_j[j] /= I;  
  
  for(j = 0; j < J; ++j){
    B += pow(x_j[j] - x__, 2);
  
    tmp = 0;
    for(i = 0; i < I; ++i)
      tmp += pow(x[i][j] - x_j[j], 2);
    tmp /= (I - 1);
    
    W += tmp;
  }
  
  B *= 1.0 / (J - 1);
  W /= J;
  
  free(x_j);
  return sqrt(1.0 - 1.0/I + (B / W) * (1.0 + 1.0 / J));
}

__host__ void oneDir(char *dir, FILE *outfp, int burnin){

  char parmName[BUF], filename[BUF], row[MAXROW];
  FILE *fp;
  int i, I, J, parmNum = 0, nRows = 0;
  num_t **x, fac;
  
  J = getNumFiles(dir);
  
  if(J < 2){
    printf("ERROR: you need at least 2 chains.\n");
    fclose(outfp);
    exit(EXIT_FAILURE);
  }
  
  sprintf(filename, "%schain0.txt", dir);
  fp = fopen(filename, "r");
  
  while(fgets(row, MAXROW, fp) != NULL)
    ++nRows;
    
  I = nRows - burnin - 1;  
  x = (num_t**) malloc(I * sizeof(num_t*));
  for(i = 0; i < I; ++i)
    x[i] = (num_t*) malloc(J * sizeof(num_t));
 
  rewind(fp);
  fscanf(fp, "%s ", parmName);

  while(!sscanf(parmName, "%f", &fac)){
    
    ++parmNum;
    
    extractOneParm(x, dir, outfp, parmName, parmNum, nRows, J, burnin); 
    fac = gelmanFactor(x, I, J);
    
    fprintf(outfp, "%s %0.3f\n", parmName, fac); 
    fscanf(fp, "%s ", parmName);
  }

  for(i = 0; i < I; ++i)
    free(x[i]);
  free(x);
    
  fclose(fp);
}

__host__ void gelmanFactors(Config *cfg){
  char file[BUF];
  FILE *fp;

  if(cfg->gelman){
  
    if(cfg->verbose)
      printf("Calculating Gelman potential scale reduction factors.\n");
  
    if(USE_R){
      system("R CMD BATCH ../R/gelman-factors.r");
      
    } else {
    
	  strcpy(file, "../out/diagnostics/gelman-factors.txt");
	  fp = fopen(file, "w");
	
	  if(fp == NULL){
		printf("ERROR: could not create Gelman factors file, \"%s\".\n", file);
		fclose(fp);
		exit(EXIT_FAILURE);
	  }
	
	  fprintf(fp, "parameter gelman-factor\n");
	
	  if(cfg->hyper)
		oneDir("../out/hyper/", fp, cfg->burnin);
	  
	  if(cfg->parms)
		oneDir("../out/parms/", fp, cfg->burnin);
	   
	  fclose(fp);
    }
  }  
}