#include <Chain.h>
#include <Config.h>
#include <constants.h>
#include <getopt.h>
#include <functions.h>
#include <stdio.h>
#include <stdlib.h>

void getopts(Config *cfg, int argc, char **argv){
  int c, option_index;
  
  struct option long_options[] = {
    {"data", required_argument, 0, 'i'},
    {"group", required_argument, 0, 'g'},
    {"probs", required_argument, 0, 'o'},  
    {"hyper", required_argument, 0, 'h'},
    {"rates", required_argument, 0, 'r'},
    {"burnin", required_argument, 0, 'b'},
    {"joint", no_argument, 0, 'j'},
    {"parms", required_argument, 0, 'p'},  
    {"sigma-c0", required_argument, 0, 'x'},
    {"d0", required_argument, 0, 'f'},
    {"a-tau", required_argument, 0, 'k'},
    {"a-alpha", required_argument, 0, 'l'},
    {"a-delta", required_argument, 0, 'm'},  
    {"b-tau", required_argument, 0, 'n'},
    {"b-alpha", required_argument, 0, 'a'},
    {"b-delta", required_argument, 0, 'c'},
    {"gamma-phi", required_argument, 0, 'q'},
    {"gamma-alpha", required_argument, 0, 'e'},  
    {"gamma-delta", required_argument, 0, 's'},
    {"sigma-phi0", required_argument, 0, 't'},
    {"sigma-alpha0", required_argument, 0, 'u'},
    {"sigma-delta0", required_argument, 0, 'v'},
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
                    "a:b:c:d:e:f:g:h:i:jk:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:1:2:3:4:5:6:7:",
                    long_options, &option_index);
    
    if(c == -1)
      break;
      
    if(c == 'i'){ /* data */
      strcpy(cfg->dataFile, optarg);
      
    } else if(c == 'g') { /* group */
      strcpy(cfg->groupFile, optarg);
            
    } else if(c == 'o') { /* probs */
      strcpy(cfg->probsFile, optarg);   
    
    } else if(c == 'h') { /* hyper */
      strcpy(cfg->hyperFile, optarg);
    
    } else if(c == 'r') { /* rates */
      strcpy(cfg->ratespFile, optarg);
    
    } else if(c == 'p') { /* parms */
      strcpy(cfg->parmsFile, optarg);
    
    } else if(c == 'b') { /* burnin */
      cfg->burnin = atoi(optarg);
    
    } else if(c == 'j') { /* joint */
      cfg->joint = atoi(optarg);
    
    } else if(c == 'x') { /* sigma-c0 */
      cfg->sigC0 = atoi(optarg);
    
    } else if(c == 'f') { /* d0 */
      cfg->d0 = atoi(optarg);
    
    } else if(c == 'k') { /* a-tau */
      cfg->aTau = atoi(optarg);
    
    } else if(c == 'l') { /* a-alpha */
      cfg->aAlp = atoi(optarg);
    
    } else if(c == 'm') { /* a-delta */
      cfg->aDel = atoi(optarg);
    
    } else if(c == 'n') { /* b-tau */
      cfg->bTau = atoi(optarg);
    
    } else if(c == 'a') { /* b-alpha */
      cfg->bAlp = atoi(optarg);
          
    } else if(c == 'c') { /* b-delta */
      cfg->bDel = atoi(optarg);    
    
    } else if(c == 'q') { /* gamma-phi */
      cfg->gamPhi = atoi(optarg);    

    } else if(c == 'e') { /* gamma-alpha */
      cfg->gamAlp = atoi(optarg);
    
    } else if(c == 's') { /* gamma-delta */
      cfg->gamDel = atoi(optarg);
    
    } else if(c == 't') { /* sigma-phi0 */
      cfg->sigPhi0 = atoi(optarg);
    
    } else if(c == 'u') { /* sigma-alpha0 */
      cfg->sigAlp0 = atoi(optarg);
    
    } else if(c == 'v') { /* sigma-delta0 */
      cfg->sigDel0 = atoi(optarg);
    
    } else if(c == 'w') { /* sigma-c */
      cfg->sigC = atoi(optarg);
      cfg->constSigC = 1;
    
    } else if(c == 'd') { /* d */
      cfg->d = atoi(optarg);
      cfg->constD = 1;
    
    } else if(c == 'y') { /* tau */
      cfg->tau = atoi(optarg);
      cfg->constTau = 1;
         
    } else if(c == 'z') { /* theta-phi */
      cfg->thePhi = atoi(optarg);
      cfg->constThePhi = 1;
    
    } else if(c == '1') { /* theta-alpha */
      cfg->theAlp = atoi(optarg);
      cfg->constTheAlp = 1;
    
    } else if(c == '2') { /* theta-delta */
      cfg->theDel = atoi(optarg);
      cfg->constTheDel = 1;
    
    } else if(c == '3') { /* sigma-phi */
      cfg->sigPhi = atoi(optarg);
      cfg->constSigPhi = 1;
    
    } else if(c == '4') { /* sigma-alpha */
      cfg->sigAlp = atoi(optarg);
      cfg->constSigAlp = 1;
    
    } else if(c == '5') { /* sigma-delta */
      cfg->sigDel = atoi(optarg);
      cfg->constSigDel = 1;
    
    } else if(c == '6') { /* pi-alpha */
      cfg->piAlp = atoi(optarg);
      cfg->constPiAlp = 1;
    
    } else if(c == '7') { /* pi-delta */
      cfg->piDel = atoi(optarg);
      cfg->constPiDel = 1;
    
    }
  }
}