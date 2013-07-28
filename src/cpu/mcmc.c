#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <functions.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>

void oneChain(Chain *a, Config *cfg){
  ++cfg->chainNum;

  if(cfg->verbose)
    printf("\n  Chain %d of %d.\n", cfg->chainNum, cfg->chains);        
        
  runChain(a, cfg);
  summarizeChain(a, cfg);
  resetChain(a, cfg);
}

void mcmc(int *pargc, char **argv){
  int i, argc = *pargc;
  Config *cfg = config(argc, argv); 
  Chain *a = newChain(cfg);

  if(cfg->debug)
    printConfig(cfg);

  if(a == NULL){
    free(cfg);
  }
  
  if(cfg->verbose)
    printf("Running %d chain(s).\n", cfg->chains);
  
  for(i = 0; i < cfg->chains; ++i)
    oneChain(a, cfg);
  
  freeChain(a, cfg);
  
  if(cfg->verbose)
    printf("Done. MCMC output written to directory: %s.\n", cfg->outDir);

  chdir(cfg->cwd);     
  free(cfg);
} 