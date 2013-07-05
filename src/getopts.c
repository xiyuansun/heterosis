#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

void getopts(int argc, char **argv){
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
      printf("data\n");
      
    } else if(c == 'g') { /* group */
      printf("group\n");
      
    } else if(c == 'o') { /* probs */
      printf("probs\n");    
    
    } else if(c == 'h') { /* hyper */
      printf("hyper\n");
    
    } else if(c == 'r') { /* rates */
      printf("rates\n");
    
    } else if(c == 'b') { /* burnin */
      printf("burnin\n");
    
    } else if(c == 'j') { /* joint */
      printf("joint\n");
    
    } else if(c == 'p') { /* parms */
      printf("parms\n");
    
    } else if(c == 'x') { /* sigma-c0 */
      printf("sigma-c0\n");
    
    } else if(c == 'f') { /* d0 */
      printf("d0\n");
    
    } else if(c == 'k') { /* a-tau */
      printf("a-tau\n");
    
    } else if(c == 'l') { /* a-alpha */
      printf("a-alpha\n");
    
    } else if(c == 'm') { /* a-delta */
      printf("a-delta\n");
    
    } else if(c == 'n') { /* b-tau */
      printf("b-tau\n");
    
    } else if(c == 'a') { /* b-alpha */
      printf("b-alpha\n");
          
    } else if(c == 'c') { /* b-delta */
      printf("b-delta\n");    
    
    } else if(c == 'q') { /* gamma-phi */
      printf("gamma-phi\n");    

    } else if(c == 'e') { /* gamma-alpha */
      printf("gamma-alpha\n");
    
    } else if(c == 's') { /* gamma-delta */
      printf("gamma-delta\n");
    
    } else if(c == 't') { /* sigma-phi0 */
      printf("sigma-phi0\n");
    
    } else if(c == 'u') { /* sigma-alpha0 */
      printf("sigma-alpha0\n");
    
    } else if(c == 'v') { /* sigma-delta0 */
      printf("sigma-delta0\n");
    
    } else if(c == 'w') { /* sigma-c */
      printf("sigma-c\n");
    
    } else if(c == 'd') { /* d */
      printf("d\n");
    
    } else if(c == 'y') { /* tau */
      printf("tau\n"); 
         
    } else if(c == 'z') { /* theta-phi */
      printf("theta-phi\n");
    
    } else if(c == '1') { /* theta-alpha */
      printf("theta-alpha\n");
    
    } else if(c == '2') { /* theta-delta */
      printf("theta-delta\n");
    
    } else if(c == '3') { /* sigma-phi */
      printf("sigma-phi\n");
    
    } else if(c == '4') { /* sigma-alpha */
      printf("sigma-alpha\n");
    
    } else if(c == '5') { /* sigma-delta */
      printf("sigma-delta\n");
    
    } else if(c == '6') { /* pi-alpha */
      printf("pi-alpha\n");
    
    } else if(c == '7') { /* pi-delta */
      printf("pi-delta\n");
    
    }
  }
}