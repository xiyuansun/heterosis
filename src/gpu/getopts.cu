#include <Config.h>
#include <constants.h>
#include <getopt.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void getopts(Config *cfg, int argc, char **argv){
  int c, option_index;
  
  struct option long_options[] = {
    {"data", required_argument, 0, 'i'},
    {"group", required_argument, 0, 'g'},
    {"chains", required_argument, 0, 'c'},
    {"iter", required_argument, 0, 'M'},
    {"burnin", required_argument, 0, 'b'},
    {"rates", no_argument, 0, 'r'},  
    {"hyper", no_argument, 0, 'h'},
    {"parms", no_argument, 0, 'p'},
    {"joint", no_argument, 0, 'j'},  
    {"seed", required_argument, 0, 's'}, 
    {"verbose", no_argument, 0, 'v'},  
    {"sigma-c0", required_argument, 0, 'x'},
    {"d0", required_argument, 0, 'f'},
    {"a-tau", required_argument, 0, 'k'},
    {"a-alpha", required_argument, 0, '0'},
    {"a-delta", required_argument, 0, '9'},  
    {"b-tau", required_argument, 0, 'n'},
    {"b-alpha", required_argument, 0, 'A'},
    {"b-delta", required_argument, 0, 'Z'},
    {"gamma-phi", required_argument, 0, 'q'},
    {"gamma-alpha", required_argument, 0, 'e'},  
    {"gamma-delta", required_argument, 0, '8'},
    {"sigma-phi0", required_argument, 0, 't'},
    {"sigma-alpha0", required_argument, 0, 'u'},
    {"sigma-delta0", required_argument, 0, 'l'},
    {"sigma-c", required_argument, 0, 'w'},
    {"d", required_argument, 0, 'd'},
    {"tau", required_argument, 0, 'y'},
    {"theta-phi", required_argument, 0, 'z'},
    {"theta-alpha", required_argument, 0, '1'},
    {"theta-delta", required_argument, 0, '2'},
    {"sigma-phi", required_argument, 0, '3'},
    {"sigma-alpha", required_argument, 0, '4'},
    {"sigma-delta", required_argument, 0, '5'}, 
    {"pi-alpha", required_argument, 0, '6'},
    {"pi-delta", required_argument, 0, '7'},                                       
    {0, 0, 0, 0}
  };

  while(1){
  
    option_index = 0;
    c = getopt_long(argc, argv, 
        "A:ab:B:c:C:d:e:f:g:G:hHi:jJk:l:m:M:n:o:pPq:rRs:S:t:u:v:w:x:y:z:1:2:3:4:5:6:7:8:9:",
        long_options, &option_index);
    
    if(c == -1)
      break;
      
    if(c == 'i' || c == 'D'){ /* data */
      strcpy(cfg->dataFile, optarg);
      
    } else if(c == 'g' || c == 'G') { /* group */
      strcpy(cfg->groupFile, optarg);
            
    } else if(c == 'c' || c == 'C') { /* number of chains */
      cfg->chains = atoi(optarg);
      
    } else if(c == 'm' || c == 'M') { /* iterations */
      cfg->M = atoi(optarg);
    
    } else if(c == 'b' || c == 'B') { /* burnin */
      cfg->burnin = atoi(optarg);
      
    } else if(c == 'r' || c == 'R') { /* rates */
      cfg->ratesFlag = 1;
    
    } else if(c == 'h' || c == 'H') { /* hyper */
      cfg->hyperFlag = 1;
    
    } else if(c == 'p' || c == 'P') { /* parms */
      cfg->parmsFlag = 1;
    
    } else if(c == 'j' || c == 'J') { /* joint */
      cfg->joint = 1;

    } else if(c == 's' || c == 'S') { /* seed */
      cfg->seed = atoi(optarg);
    
    } else if(c == 'v' || c == 'V') { /* verbose */
      cfg->verbose = 1;
    
    } else if(c == 'x') { /* sigma-c0 */
      cfg->sigC0 = atof(optarg);
    
    } else if(c == 'f') { /* d0 */
      cfg->d0 = atof(optarg);
    
    } else if(c == 'k') { /* a-tau */
      cfg->aTau = atof(optarg);
    
    } else if(c == '0') { /* a-alpha */
      cfg->aAlp = atof(optarg);
    
    } else if(c == '9') { /* a-delta */
      cfg->aDel = atof(optarg);
    
    } else if(c == 'n') { /* b-tau */
      cfg->bTau = atof(optarg);
    
    } else if(c == 'A') { /* b-alpha */
      cfg->bAlp = atof(optarg);
          
    } else if(c == 'Z') { /* b-delta */
      cfg->bDel = atof(optarg);    
    
    } else if(c == 'q') { /* gamma-phi */
      cfg->gamPhi = atof(optarg);    

    } else if(c == 'e') { /* gamma-alpha */
      cfg->gamAlp = atof(optarg);
    
    } else if(c == '8') { /* gamma-delta */
      cfg->gamDel = atof(optarg);
    
    } else if(c == 't') { /* sigma-phi0 */
      cfg->sigPhi0 = atof(optarg);
    
    } else if(c == 'u') { /* sigma-alpha0 */
      cfg->sigAlp0 = atof(optarg);
    
    } else if(c == 'l') { /* sigma-delta0 */
      cfg->sigDel0 = atof(optarg);
    
    } else if(c == 'w') { /* sigma-c */
      cfg->sigC = atof(optarg);
      cfg->constSigC = 1;
    
    } else if(c == 'd') { /* d */
      cfg->d = atof(optarg);
      cfg->constD = 1;
    
    } else if(c == 'y') { /* tau */
      cfg->tau = atof(optarg);
      cfg->constTau = 1;
         
    } else if(c == 'z') { /* theta-phi */
      cfg->thePhi = atof(optarg);
      cfg->constThePhi = 1;
    
    } else if(c == '1') { /* theta-alpha */
      cfg->theAlp = atof(optarg);
      cfg->constTheAlp = 1;
    
    } else if(c == '2') { /* theta-delta */
      cfg->theDel = atof(optarg);
      cfg->constTheDel = 1;
    
    } else if(c == '3') { /* sigma-phi */
      cfg->sigPhi = atof(optarg);
      cfg->constSigPhi = 1;
    
    } else if(c == '4') { /* sigma-alpha */
      cfg->sigAlp = atof(optarg);
      cfg->constSigAlp = 1;
    
    } else if(c == '5') { /* sigma-delta */
      cfg->sigDel = atof(optarg);
      cfg->constSigDel = 1;
    
    } else if(c == '6') { /* pi-alpha */
      cfg->piAlp = atof(optarg);
      cfg->constPiAlp = 1;
    
    } else if(c == '7') { /* pi-delta */
      cfg->piDel = atof(optarg);
      cfg->constPiDel = 1;
    
    } else { /* error */
      printf("Argument error. See usage.\n");
      free(cfg);
      exit(EXIT_FAILURE);
    
    }
  }
}