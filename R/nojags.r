# Author: Will Landau
# PIs: Dan Nettleton, Jarad Niemi, Peng Liu
# Iowa State University
# May 2013
# 
# This code attempts to learn about differential expression
# between 2 treatment groups of an RNA-Seq dataset. See
# writeup.pdf for the model.

library(Biobase)
library(coda)

# sampling from known distributions

sampleNormal = function(m = 0, s = 1){ # device
  u1 = runif(1);
  u2 = runif(1);
  return(sqrt(-2 * log(u1)) * sin(2 * 3.14159265 * u2) * s + m);
}

sampleGamma = function(shape = 1, rate = 1, lb = 0){ # device
  if(shape <= 0){
    print(paste("Error: bad shape:", shape));
    return(0/0);
  }
  
  if(rate <= 0){
    print(paste("Error: bad rate:", rate));
    return(0/0);
  }

  if(shape - 1 < lb * rate){ # Chung (1998)
    c = lb * rate;

    eps0 = (c - shape + sqrt((c - shape)^2 + 4 * c))/(2 * c);

    if(eps0 > 1){
      eps = 0.75;
    } else {
      eps = eps0
    }

    A = ((1 - eps)/(shape - 1))^(shape - 1) * exp(shape - 1);

    while(1){
      x = - (1/eps) * log(runif(1)) + lb * rate;
      u = runif(1);

      if(u < A * x^(shape - 1) * exp((eps - 1) * x))
        return(x / rate);

    }

  } else if(shape >= 1){ # Marsaglia and Tsang (2000)

    d = shape - 1/3;
    c = 1 / sqrt(9 * d);

    while(1){
      v = -1;
      while(v <= 0){
        x = sampleNormal(0, 1);
        v = (1 + c*x)^3;
      }

      ret = d * v / rate

      if(ret > lb){
        u = runif(1);

        if(u < 1 - 0.0331 * x^4)
          return(ret);

        if(log(u) < 0.5 * x^2 + d * (1 - v + log(v)))
          return(ret);
      }
    }
  } else if (0.135 <= shape && shape < 1){ # Kundu and Gupta (2006)

    while(1){      

      u = runif(1);
      x = -2 * log(1 - u^(1 / shape));
      ret = x / rate;

      if(ret > lb){
        v = runif(1);

        tmp1 = exp(-x/2);
        tmp2 = x^(shape - 1)* tmp1 * 2^(1 - shape) * (1 - tmp1)^(1 - shape);

        if(v < tmp2)
          return(ret);
      }
    }
  } else{ # Martin and Liu (2013)
   
    while(1){ # 
      lam = 1/shape - 1;
      w = shape / (exp(1 - shape));
      r = 1/(1 + w); 

      u = runif(1);

      if(u <= r){
        z = - log(u/r);
      } else {
        z = log(runif(1))/lam;
      }
      
      ret = exp(-z/shape) / rate;

      if(ret > lb){
        if(z >= 0){
          haznaz = exp(-exp(-z / shape));
        } else{
          haznaz = 1/(w * lam) * exp((lam - 1) * z -exp(z / shape));
        }

        if(haznaz > runif(1))
          return(ret);
      }
    }
  }
}

sampleBeta = function(a, b){ # device
  x = sampleGamma(a, 1, 0);
  y = sampleGamma(b, 1, 0);
  return(x / (x + y));
}

# data structures and initialization

allocConfig = function(){
  list(
    datafile = "",
    groupfile = "",

    heterosisfile = "",
    diffexprfile = "",

    hyperfile = "",
    ratesfile = "",

    iter = 0,
    burnin = 0,
    joint = 0,

    sigC0 = 0,
    d0 = 0,
    
    aTau = 0,
    aAlp = 0,
    aDel = 0,
    
    bTau = 0,
    bAlp = 0,
    bDel = 0,
    
    gamPhi = 0,
    gamAlp = 0,
    gamDel = 0,
    
    sigPhi0 = 0,
    sigAlp0 = 0,
    sigDel0 = 0,
    
    sigC = 0,
    
    d = 0,
    tau = 0,
    
    thePhi = 0,
    theAlp = 0,
    theDel = 0,

    sigPhi = 0,
    sigAlp = 0,
    sigDel = 0,

    piAlp = 0,
    piDel = 0,

    skip = list(

      sigC = 0,
    
      d = 0,
      tau = 0,
    
      thePhi = 0,
      theAlp = 0,
      theDel = 0,

      sigPhi = 0,
      sigAlp = 0,
      sigDel = 0,

      piAlp = 0,
      piDel = 0, 
      
    )
  );
}

getopts = function(){
  cfg = allcoConfig();

  # initialization constants = -1 when unset

  cfg$sigC0 = -1; 
  cfg$d0 = -1; 

  cfg$aTau = -1; 
  cfg$aAlp = -1; 
  cfg$aDel = -1; 

  cfg$bTau = -1; 
  cfg$bAlp = -1; 
  cfg$bDel = -1; 
  
  cfg$gamPhi = -1; 
  cfg$gamAlp = -1; 
  cfg$gamDel = -1; 

  cfg$sigPhi0 = -1; 
  cfg$sigAlp0 = -1; 
  cfg$sigDel0 = -1; 

  cfg$skip$sigC = 0;

  cfg$skip$d = 0;
  cfg$skip$tau = 0;

  cfg$skip$thePhi = 0;
  cfg$skip$theAlp = 0;
  cfg$skip$theDel = 0; 

  cfg$skip$sigPhi = 0; 
  cfg$skip$sigAlp = 0; 
  cfg$skip$sigDel = 0; 

  cfg$skip$piAlp = 0;
  cfg$skip$piDel = 0;

  # Use long getopt to parse command line args.

  cfg$datafile = "~/heterosis/var/data/hammer/hammer.txt";
  cfg$groupfile = "~/heterosis/var/data/hammer/group.txt"; 

  heterosisfile = "";
  diffexprfile = "";

  hyperfile = "";
  ratesfile = "";

  iter = 50;
  burnin = 0;
  joint = 0;

  if(cfg$sigC0 < 0)
    cfg$sigC0 = 10; 

  if(cfg$d0 < 0)
    cfg$d0 = 1e3; 

  if(cfg$aTau < 0)
    cfg$aTau = 1e2; 

  if(cfg$aAlp < 0)
    cfg$aAlp = 1; 

  if(cfg$aDel < 0)
    cfg$aDel = 1; 

  if(cfg$bTau < 0)
    cfg$bTau = 1e2; 

  if(cfg$bAlp < 0)
    cfg$bAlp = 1; 

  if(cfg$bDel < 0)
    cfg$bDel = 1; 
  
  if(cfg$gamPhi < 0)
    cfg$gamPhi = 2; 

  if(cfg$gamAlp < 0)
    cfg$gamAlp = 2; 

  if(cfg$gamDel < 0)
    cfg$gamDel = 2; 

  if(cfg$sigPhi0 < 0)
    cfg$sigPhi0 = 2; 

  if(cfg$sigAlp0 < 0)
    cfg$sigAlp0 = 1e2;

  if(cfg$sigDel0 < 0)
    cfg$sigDel0 = 1e2; 

  if(!cfg$skip$sigC)
    cfg$sigC = runif(1, 0, cfg$sigC0);

  if(!cfg$skip$d)
    cfg$d = runif(1, 0, cfg$d0);

  if(!cfg$skip$tau)
    cfg$tau = sqrt(sampleGamma(shape = cfg$aTau, rate = cfg$bTau));

  if(!cfg$skip$thePhi)
    cfg$thePhi = sampleNormal(0, cfg$gamPhi);

  if(!cfg$skip$theAlp)
    cfg$theAlp = sampleNormal(0, cfg$gamAlp);

  if(!cfg$skip$theDel)
    cfg$theDel = sampleNormal(0, cfg$gamDel); 

  if(!cfg$skip$sigPhi)
    cfg$sigPhi = runif(1, 0, cfg$sigPhi0); 

  if(!cfg$skip$sigAlp)
    cfg$sigAlp = runif(1, 0, cfg$sigAlp0); 

  if(!cfg$skip$sigDel)
    cfg$sigDel = runif(1, 0, cfg$sigDel0); 

  if(!cfg$skip$piAlp)
    cfg$piAlp = sampleBeta(cfg$aAlp, cfg$bAlp);

  if(!cfg$skip$piDel)
    cfg$piDel = sampleBeta(cfg$aDel, cfg$bDel);
 
  cfg
}

allocChain = function(M, N, G){ # host (bunch of cudaMallocs)
  list(

    # data 

    y = array(0, c(N, G)),
    yMeanG = rep(0, N),
    grp = rep(0, N),

    M = M,
    N = N,
    G = G,
    
    burnin = 0, # burn-in for calculating heterosis 
                # and differential expression probabilities

    # initialization constants

    sigC0 = 0,
    d0 = 0,

    aTau = 0,
    aAlp = 0,
    aDel = 0,

    bTau = 0,
    bAlp = 0,
    bDel = 0,
  
    gamPhi = 0,
    gamAlp = 0,
    gamDel = 0,

    sigPhi0 = 0,
    sigAlp0 = 0,
    sigDel0 = 0,

    # parameters

    c = array(0, c(M + 1, N)),
      sigC = rep(0, M + 1),
    
    eps = array(0, c(M + 1, N, G)),
      eta = array(0, c(M + 1, G)),
        d = rep(0, M + 1),
        tau = rep(0, M + 1),

    phi = array(0, c(M + 1, G)),
      thePhi = rep(0, M + 1),
      sigPhi = rep(0, M + 1),

    alp = array(0, c(M + 1, G)),
      theAlp = rep(0, M + 1),
      sigAlp = rep(0, M + 1),
      piAlp = rep(0, M + 1),

    del = array(0, c(M + 1, G)),
      theDel = rep(0, M + 1),
      sigDel = rep(0, M + 1),
      piDel = rep(0, M + 1),

    # temporary and return values

    tmp1 = rep(0, G),
    tmp2 = rep(0, G),  

    shape = 0,
    rate = 0,
    
    old = array(0, c(N, G)),
    new = array(0, c(N, G)),

    lOld = array(0, c(N, G)),
    lNew = array(0, c(N, G)),

    # struct to store the current place in the chain for each parameter

    m = list(
      c = 0,
      sigC = 0,

      eps = 0,
      eta = 0,
      d = 0,
      tau = 0,

      phi = 0,
      alp = 0,
      del = 0,

      thePhi = 0,
      theAlp = 0,
      theDel = 0,

      sigPhi = 0,
      sigAlp = 0,
      sigDel = 0,

      piAlp = 0,
      piDel = 0
    ),

    # tuning parameters for metropolis steps (std deviations of normal distributions)

    tun = list(
      c = rep(0, N),
      eps = array(0, c(N, G)),
      d = 0,
      phi = rep(0, G)
    ),

    # number of acceptances for metropolis steps

    acc = list(
      c = rep(0, N),
      eps = array(0, c(N, G)),
      d = 0,
      phi = rep(0, G),
      alp = rep(0, G),
      del = rep(0, G)
    )
  );
}

newChain_kernel1 = function(a){ # kernel: G blocks, 1 thread each
 for(g in 1:a$G){
    a$phi[1, g] = sampleNormal(a$thePhi[1], a$sigPhi[1]);

    u = runif(1);
    if(u < a$piAlp[1]){
      a$alp[1, g] = 0;
    } else {
      a$alp[1, g] = sampleNormal(a$theAlp[1], a$sigAlp[1]);
    }
    
    u = runif(1);
    if(u < a$piDel[1]){
      a$del[1, g] = 0;
    } else {
      a$del[1, g] = sampleNormal(a$theDel[1], a$sigDel[1]);
    }
 
    a$eta[1, g] = 1/sqrt(sampleGamma(shape = a$d[1] / 2, 
                                  rate = a$d[1] * a$tau[1]^2 / 2))

    for(n in 1:a$N)
      a$eps[1, n, g] = sampleNormal(0, a$eta[1, g]);
  }

  return(a);
}

newChain_kernel2 = function(a){ # kernel: 1 block, 1 thread

  # location in chain:

  a$m$c = 1;
  a$m$sigC = 1;

  a$m$eps = 1;
  a$m$eta = 1;
  a$m$d = 1;
  a$m$tau = 1;

  a$m$phi = 1;
  a$m$alp = 1;
  a$m$del = 1;

  a$m$thePhi = 1;
  a$m$theAlp = 1;
  a$m$theDel = 1;

  a$m$sigPhi = 1;
  a$m$sigAlp = 1;
  a$m$sigDel = 1;

  a$m$piAlp = 1;
  a$m$piDel = 1;

  # tuning parameters for metropolis steps (std deviations of normal dists): 

  a$tun$d = 500;

  for(n in 1:a$N)
    a$tun$c[n] = 1;

  for(g in 1:a$G){
    a$tun$phi[g] = 1;

    for(n in 1:a$N)
      a$tun$eps[n, g] = 1;
  }

  # number of acceptances for metropolis steps

  a$acc$d = 0;

  for(n in 1:a$N){
    a$acc$c[n] = 0;
  
    for(g in 1:a$G)
      a$acc$eps[n, g] = 0;
  }

  for(g in 1:a$G){
    a$acc$phi[g] = 0;
    a$acc$alp[g] = 0;
    a$acc$del[g] = 0;
  }
  
  a
}

newChain = function(cfg){ # host (bunch of cudaMemcpies and kernels)

  y = t(matrix(read.table(cfg$datafile)))
  group = scan(cfg$groupfile)

  M = cfg$iter;
  N = dim(y)[1];
  G = dim(y)[2];

  a = allocChain(M, N, G);
  
  # data

  for(n in 1:N){
    a$grp[n] = grp[n]; # CudaMemcpy
    tmp = 0;

    for(g in 1:G){
      a$y[n, g] = y[n, g]; # CudaMemcpy
      tmp = tmp + y[n, g];
    }

    a$yMeanG[n] = tmp / G; # CudaMemcpy
  }

  a$M = M; # CudaMemcpy
  a$N = N; # CudaMemcpy
  a$G = G; # CudaMemcpy

  # burn-in for calculating heterosis 
  # and differential expression probabilities

  a$burnin = M/2; 

  # initialization constants: CudaMemcpies

  a$sigC0 = cfg$sigC0; # CudaMemcpy
  a$d0 = cfg$d0; # CudaMemcpy

  a$aTau = cfg$aTau; # CudaMemcpy
  a$aAlp = cfg$aAlp; # CudaMemcpy
  a$aDel = cfg$aDel; # CudaMemcpy

  a$bTau = cfg$bTau; # CudaMemcpy
  a$bAlp = cfg$bAlp; # CudaMemcpy
  a$bDel = cfg$bDel; # CudaMemcpy
  
  a$gamPhi = cfg$gamPhi; # CudaMemcpy
  a$gamAlp = cfg$gamAlp; # CudaMemcpy
  a$gamDel = cfg$gamDel; # CudaMemcpy

  a$sigPhi0 = cfg$sigPhi0; # CudaMemcpy
  a$sigAlp0 = cfg$sigAlp0; # CudaMemcpy
  a$sigDel0 = cfg$sigDel0; # CudaMemcpy

  # hyperparameters: CudaMemcpies

  a$sigC[1] = cfg$sigC # CudaMemcpy

  a$d[1] = cfg$d # CudaMemcpy
  a$tau[1] = cfg$tau # CudaMemcpy

  a$thePhi[1] = cfg$thePhi # CudaMemcpy
  a$theAlp[1] = cfg$theAlp # CudaMemcpy
  a$theDel[1] = cfg$theDel # CudaMemcpy

  a$sigPhi[1] = cfg$sigPhi # CudaMemcpy
  a$sigAlp[1] = cfg$sigAlp # CudaMemcpy
  a$sigDel[1] = cfg$sigDel # CudaMemcpy

  a$piAlp[1] = cfg$piAlp # CudaMemcpy
  a$piDel[1] = cfg$piDel # CudaMemcpy

  # compute initial normalization factors, mostly using priors

  lqts = rep(0, N);
  tmp = rep(0, G);

  s = 0;
  for(n in 1:N){
    for(g in 1:G)
      tmp[g] = y[n, g];

    tmp = sort(tmp); # use qsort()
    lqts[n] = log(tmp[floor(G * 0.75)]); 
    s = s + lqts[n];
  }

  s = s / N;

  for(n in 1:N)
    a$c[1, n] = lqts[n] - s; # CudaMemcpy

  a = newChain_kernel1(a);
  a = newChain_kernel2(a);

  return(a);
}

mu = function(a, n, phi, alp, del){ # device
  if(a$grp[n] == 1){
    return(phi - alp);
  } else if(a$grp[n] == 2){
    return(phi + del);
  } else if(a$grp[n] == 3){
    return(phi + alp);
  }
}


# log full conditionals with no convenient form

lC_kernel1 = function(a, n){ # kernel <<<G, 1>>>
  for(g in 1:a$G) # PARALLELIZE
    a$tmp1[g] = exp(a$eps[a$m$eps, n, g] + mu(a, n, a$phi[a$m$phi, g], 
                                            a$alp[a$m$alp, g], a$del[a$m$del, g]));

  a
}

lC_kernel2 = function(a, n){ # kernel <<<G, 1>>>
  a$tmp2[1] = 0;
  for(g in 1:a$G) # PARALLEL PAIRWISE SUM IN THRUST
    a$tmp2[1] = a$tmp2[1] + a$tmp1[g];

  a
}

lC_kernel3 = function(a, n, newArg){ # kernel <<<1, 1>>>

  if(newArg){
    arg = a$new[1, 1];
  } else {
    arg = a$old[1, 1];
  }

  ret = arg * a$G * a$yMeanG[n] - exp(arg) * a$tmp2[1] - (arg*arg) / 
        (2 * a$sigC[a$m$sigC] * a$sigC[a$m$sigC]);

  if(newArg){
    a$lNew[1, 1] = ret;
  } else {
    a$lOld[1, 1] = ret;
  }

  a
}

lC = function(a, n, newArg){ # host
  a = lC_kernel1(a, n);
  a = lC_kernel2(a, n);
  a = lC_kernel3(a, n, newArg);
  a
}

lEps = function(a, n, g, arg){ # device

  ret = a$y[n, g] * arg - 
        exp(a$c[a$m$c, n] + arg + mu(a, n, a$phi[a$m$phi, g], 
                                     a$alp[a$m$alp, g], a$del[a$m$del, g])) 
          - arg^2 / (2 * a$eta[a$m$eta, g]^2);

  return(ret);
}

lD_kernel1 = function(a){ # kernel <<<G, 1>>>
  for(g in 1:a$G){ 
    a$tmp1[g] = 2 * log(a$eta[a$m$eta, g]);
    a$tmp2[g] = 1/(a$eta[a$m$eta, g] * a$eta[a$m$eta, g]);
  }
  a
}

lD_kernel2 = function(a){ # kernel: pairwise sum in Thrust
  a$s1 = 0;
  for(g in 1:a$G) # PARALLELIZE
    a$s1 = a$s1 + a$tmp1[g];

  a
}

lD_kernel3 = function(a){ # kernel: pairwise sum in Thrust
  a$s2 = 0;
  for(g in 1:a$G) # PARALLELIZE
    a$s2 = a$s2 + a$tmp2[g];
  a
}

lD_kernel4 = function(a, newArg){ # kernel <<<1, 1>>>
 
  if(newArg){
    arg = a$new[1, 1];
  } else{
    arg = a$old[1, 1];
  }

  a$tmp1[1] = arg * a$tau[a$m$tau]^2 / 2;

  ret = -a$G * lgamma(arg/2) + (a$G * arg / 2) * log(a$tmp1[1]);
  ret = ret  - (arg/2 + 1) * a$s1 - a$tmp1[1] * a$s2;

  if(newArg){
    a$lNew[1, 1] = ret;
  } else{
    a$lOld[1, 1] = ret;
  }

  a
}

lD = function(a, newArg){ # host
  
  # use cudamemcpies
  if(newArg){
    if(a$new[1, 1] <= 0 || a$new[1, 1] > a$d0){
      a$lNew[1, 1] = -Inf;
      return(a);
    }
  } else {
    if(a$old[1, 1] <= 0 || a$old[1, 1] > a$d0){
      a$lOld[1, 1] = -Inf; 
      return(a);
    }
  }

  a = lD_kernel1(a);
  a = lD_kernel2(a);
  a = lD_kernel3(a);
  a = lD_kernel4(a, newArg)

  a
}

lPhi = function(a, g, arg){ # device
 
  s = 0; 
  tmp = 0;

  for(n in 1:a$N){
    tmp = mu(a, n, arg, a$alp[a$m$alp, g], a$del[a$m$del, g]);
    s = s + a$y[n, g] * tmp - exp(a$c[a$m$c, n] + 
        a$eps[a$m$eps, n, g] + tmp);
  }
 
  ret = s - (arg - a$thePhi[a$m$thePhi])^2 / (2 * a$sigPhi[a$m$sigPhi]^2);
  return(ret);
}

lAlp = function(a, g, arg){ # device
  
  s = 0; 
  for(n in 1:a$N){
    if(a$grp[n] != 2){
      tmp = mu(a, n, a$phi[a$m$phi, g], arg, a$del[a$m$del, g]);
      s = s + a$y[n, g] * tmp - exp(a$c[a$m$c, n] + 
          a$eps[a$m$eps, n, g] + tmp);
    }
  }
 
  if(arg * arg > 1e-6){
    tmp = -(arg - a$theAlp[a$m$theAlp])^2 / (2 * a$sigAlp[a$m$sigAlp]^2) -
                log(1 - a$piAlp[a$m$piAlp]);
  } else {
    tmp = log(a$piAlp[a$m$piAlp]);
  }

  ret = s + tmp;
  return(ret);
}

lDel = function(a, g, arg){ # device 
  
  s = 0; 
  for(n in 1:a$N){
    if(a$grp[n] != 2){
      tmp = mu(a, n, a$phi[a$m$phi, g], a$alp[a$m$alp, g], arg);
      s = s + a$y[n, g] * tmp - exp(a$c[a$m$c, n] + 
          a$eps[a$m$eps, n, g] + tmp);
    }
  }
 
  if(arg * arg > 1e-6){
    tmp = -(arg - a$theDel[a$m$theDel])^2 / (2 * a$sigDel[a$m$sigDel]^2) -
                log(1 - a$piDel[a$m$piDel]);
  } else {
    tmp = log(a$piDel[a$m$piDel]);
  }

  ret = s + tmp;
  return(ret);
}

lPhiAlpDel = function(a, g, argPhi, argAlp, argDel){ # device
 
  s = 0; 
  tmp = 0;

  for(n in 1:a$N){
    tmp = mu(a, n, argPhi, argAlp, argDel);
    s = s + a$y[n, g] * tmp - exp(a$c[a$m$c, n] + 
        a$eps[a$m$eps, n, g] + tmp);
  }

  # phi part 
  ret = s - (argPhi - a$thePhi[a$m$thePhi])^2 / (2 * a$sigPhi[a$m$sigPhi]^2);

  # alpha part
  if(argAlp * argAlp > 1e-6){
    tmp = -(argAlp - a$theAlp[a$m$theAlp])^2 / (2 * a$sigAlp[a$m$sigAlp]^2) -
                log(1 - a$piAlp[a$m$piAlp]);
  } else {
    tmp = log(a$piAlp[a$m$piAlp]);
  }

  ret = ret + tmp;

  # delta part
  if(argDel * argDel > 1e-6){
    tmp = -(argDel - a$theDel[a$m$theDel])^2 / (2 * a$sigDel[a$m$sigDel]^2) -
                log(1 - a$piDel[a$m$piDel]);
  } else {
    tmp = log(a$piDel[a$m$piDel]);
  }

  ret = ret + tmp

  return(ret);
}



# samplers

sampleC_kernel1 = function(a, n){ # kernel <<<1, 1>>>
  a$old[1, 1] = a$c[a$m$c, n];
  a$new[1, 1] = sampleNormal(a$old[1, 1], a$tun$c[n]);

  a
}

sampleC_kernel2 = function(a, n){ # kernel <<<1, 1>>>
  lp = min(0, a$lNew[1, 1] - a$lOld[1, 1]);
  lu = log(runif(1));
    
  if(lu < lp){ # accept
    a$c[a$m$c + 1, n] = a$new[1, 1];
    a$tun$c[n] = a$tun$c[n] * 1.1; # Increase the proposal variance to avoid  
                                   # getting stuck in a mode
    a$acc$c[n] = a$acc$c[n] + 1;
  } else { # reject
    a$c[a$m$c + 1, n] = a$old[1, 1];
    a$tun$c[n] = a$tun$c[n] / 1.1; # If you're rejecting too often, decrease the  
                                   # proposal variance to sample closer to 
                                   # the last accepted value.
  }

  a
}

sampleC_kernel3 = function(a){ # kernel <<<1, 1>>>
  a$m$c = a$m$c + 1;
  a
}

sampleC = function(a){ # host
  for(n in 1:a$N){ 

    a = sampleC_kernel1(a, n);

    a = lC(a, n, newArg = 1);
    a = lC(a, n, newArg = 0);

    a = sampleC_kernel2(a, n);
  }

  a = sampleC_kernel3(a);
  a
}

sampleSigC = function(a){ # kernel <<<1, 1>>>

  for(n in 1:a$N) 
    a$tmp1[n] = a$c[a$m$c, n]^2;

  rate = 0;
  for(n in 1:a$N) 
    rate = rate + a$tmp1[n];
  
  shape = (a$N - 1) / 2; 
  rate = rate / 2;
  lb = 1 / a$sigC0^2  

  if(shape >= 1 && rate > 0){
    a$sigC[a$m$sigC + 1] = 1/sqrt(sampleGamma(shape, rate, lb));
  } else {
    a$sigC[a$m$sigC + 1] = a$sigC[a$m$sigC];
  }

  a$m$sigC = a$m$sigC + 1;
  return(a)
}

sampleEps_kernel1 = function(a){ # kernel <<<N, G>>>
  for(g in 1:a$G){
    for(n in 1:a$N){ 
      old = a$eps[a$m$eps, n, g];
      new = sampleNormal(old, a$tun$eps[n, g]);

      lp = min(0, lEps(a, n, g, new) - lEps(a, n, g, old));
      lu = log(runif(1));
      
      if(lu < lp){ # accept
        a$eps[a$m$eps + 1, n, g] = new;
        a$tun$eps[n, g] = a$tun$eps[n, g] * 1.1;
        a$acc$eps[n, g] = a$acc$eps[n, g] + 1; 
      } else { # reject
        a$eps[a$m$eps + 1, n, g] = old;
        a$tun$eps[n, g] = a$tun$eps[n, g] / 1.1;
      }
    }
  }

  a
}

sampleEps_kernel2 = function(a){ # kernel <<<1, 1>>>
  a$m$eps = a$m$eps + 1;
  a
}

sampleEps = function(a){ # host

  a = sampleEps_kernel1(a);
  a = sampleEps_kernel2(a);

  a
}

sampleEta_kernel1 = function(a){ # kernel <<<1, 1>>>
  a$shape = (a$N + a$d[a$m$d]) / 2; 
  a
}

sampleEta_kernel2 = function(a){ # kernel <<<G, 1>>>
  for(g in 1:a$G){

    for(n in 1:a$N) 
      a$tmp1[n] = a$eps[a$m$eps, n, g]^2;

    rate = 0;
    for(n in 1:a$N) 
      rate = rate + a$tmp1[n];
  
    rate = (rate + a$d[a$m$d] * a$tau[a$m$tau] * a$tau[a$m$tau]) / 2; 

    if(a$shape >= 1 && rate > 0){
      a$eta[a$m$eta + 1, g] = 1/sqrt(sampleGamma(a$shape, rate));
    } else {
      a$eta[a$m$eta + 1, g] = a$eta[a$m$eta, g];
    }
  }
  a
}

sampleEta_kernel3 = function(a){ # kernel <<<1, 1>>>
  a$m$eta = a$m$eta + 1;
  a
}

sampleEta = function(a){

  a = sampleEta_kernel1(a);
  a = sampleEta_kernel2(a);
  a = sampleEta_kernel3(a);

  return(a)
}


sampleD_kernel1 = function(a){ # kernel <<<1, 1>>>
  a$old[1, 1] = a$d[a$m$d];
  a$new[1, 1] = sampleNormal(a$old[1, 1], a$tun$d);
  a
}

sampleD_kernel2 = function(a){ # kernel <<<1, 1>>>
  lp = min(0, a$lNew[1, 1] - a$lOld[1, 1]);
  lu = log(runif(1));

  if(lu < lp){ # accept
    a$d[a$m$d + 1] = a$new[1, 1];
    a$tun$d = a$tun$d * 1.1; # Increase the proposal variance to avoid getting 
                                 # stuck in a mode
    a$acc$d = a$acc$d + 1;
  } else { # reject
    a$d[a$m$d + 1] = a$old[1, 1];
    a$tun$d = a$tun$d / 1.1; # If you're rejecting too often, decrease the proposal 
                                 # variance to sample closer to the last accepted value.
  }

  a$m$d = a$m$d + 1;
  a
}

sampleD = function(a){ 

  a = sampleD_kernel1(a);

  a = lD(a, 1);
  a = lD(a, 0);

  a = sampleD_kernel2(a);
  a
}

sampleTau_kernel1 = function(a){ # kernel <<<G, 1>>>
  for(g in 1:a$G) # PARALLELIZE
    a$tmp1[g] = 1/a$eta[a$m$eta, g]^2;

  a
}

sampleTau_kernel2 = function(a){ # pairwise sum in Thrust
  tmp = 0;
  for(g in 1:a$G) 
    tmp = tmp + a$tmp1[g];

  a$rate = tmp;
  a
}

sampleTau_kernel3 = function(a){
  rate = a$rate * a$d[a$m$d] / 2 + a$bTau;
  shape = a$aTau + a$G * a$d[a$m$d] / 2;

  if(shape >= 1 && rate > 0){
    a$tau[a$m$tau + 1] = 1/sqrt(sampleGamma(shape, rate));
  } else {
    a$tau[a$m$tau + 1] = a$tau[a$m$tau];
  }

  a$m$tau = a$m$tau + 1;

  a
}

sampleTau = function(a){ # host
  a = sampleTau_kernel1(a);
  a = sampleTau_kernel2(a);
  a = sampleTau_kernel3(a);
 
  return(a)
}

samplePhi_kernel1 = function(a){ # kernel <<<G, 1>>>
  for(g in 1:a$G){ # PARALLELIZE

    old = a$phi[a$m$phi, g];
    new = sampleNormal(old, a$tun$phi[g]);

    lp = min(0, lPhi(a, g, new) - lPhi(a, g, old));
    lu = log(runif(1));
    
    if(lu < lp){ # accept
      a$phi[a$m$phi + 1, g] = new;
      a$tun$phi[g] = a$tun$phi[g] * 1.1; 
      a$acc$phi[g] = a$acc$phi[g] + 1;
    } else { # reject
      a$phi[a$m$phi + 1, g] = old;
      a$tun$phi[g] = a$tun$phi[g] / 1.1; 
    }
  }

  a
}

samplePhi_kernel2 = function(a){ # kernel <<<1, 1>>>
  a$m$phi = a$m$phi + 1;
  a
}

samplePhi = function(a){ # host
  a = samplePhi_kernel1(a);
  a = samplePhi_kernel2(a);
  return(a)
}

sampleThePhi_kernel1 = function(a){ # pairwise sum in Thrust
  a$s1 = 0; 
  for(g in 1:a$G) # PARALLELIZE
    a$s1 = a$s1 + a$phi[a$m$phi, g];

  a
}

sampleThePhi_kernel2 = function(a){ # kernel <<<1, 1>>>
  gs = a$gamPhi^2;
  ss = a$sigPhi[a$m$sigPhi]^2;
  den = (a$G * gs + ss);

  m = gs * a$s1 / den;
  s = gs * ss / den;

  a$thePhi[a$m$thePhi + 1] = sampleNormal(m, s);
  a$m$thePhi = a$m$thePhi + 1;

  a
}

sampleThePhi = function(a){ # host
  a = sampleThePhi_kernel1(a);
  a = sampleThePhi_kernel2(a);
  return(a)
}

sampleSigPhi_kernel1 = function(a){ # kernel <<<G, 1>>>
  for(g in 1:a$G) # PARALLELIZE
    a$tmp1[g] = (a$phi[a$m$phi, g] - a$thePhi[a$m$thePhi])^2;

  a
}

sampleSigPhi_kernel2 = function(a){ # parallel pairwise sum in Thrust
  rate = 0;
  for(g in 1:a$G) # PARALLELIZE
    rate = rate + a$tmp1[g];
  a$rate = rate;  

  a
}

sampleSigPhi_kernel3 = function(a){ # kernel <<<1, 1>>>
  rate = a$rate / 2;
  shape = (a$G - 1) / 2;
  lb = 1/a$sigPhi0^2;

  if(shape >= 1 && rate > 0){
    a$sigPhi[a$m$sigPhi + 1] = 1/sqrt(sampleGamma(shape, rate, lb));
  } else {
    a$sigPhi[a$m$sigPhi + 1] = a$sigPhi[a$m$sigPhi];
  }

  a$m$sigPhi = a$m$sigPhi + 1;

  a
}

sampleSigPhi = function(a){ # host
 
  a = sampleSigPhi_kernel1(a);
  a = sampleSigPhi_kernel2(a);
  a = sampleSigPhi_kernel3(a);

  return(a)
}


piAlpPrime = function(a, g, avg, s){ # device
  prod = 1;

  for(n in 1:a$N){
    if(a$grp[n] != 2){
      tmp = a$c[a$m$c, n] + a$eps[a$m$eps, n, g] +
            mu(a, n, a$phi[a$m$phi, g], 0, a$del[a$m$del, g]);

      den = exp(a$y[n, g] * tmp - exp(tmp));

      i = a$c[a$m$c, n] + a$eps[a$m$eps, n, g] + a$phi[a$m$phi, g]

      A = -1/(2 * s * s) - 1/2;

      if(a$grp[n] == 1){
        B = i + avg/(s * s) - a$y[n, g] + 1;
      } else if(a$grp[n] == 3){
        B = -i + avg/(s * s) + a$y[n, g] - 1;
      }     

      C = -(i * i)/2 + i * a$y[n, g] - i - (avg * avg)/(2 * s * s) - 1;
      D = 1.0/sqrt(2 * pi * s * s);
                     
      num = D * exp(C - (B * B)/(4 * A)) * sqrt(- pi/A);

      if(den != 0 && is.finite(num))
        prod = prod * num/den

    }
  }

  ret = ((1 - a$piAlp[a$m$piAlp]) / a$piAlp[a$m$piAlp]) * prod;
  ret = 1/(1 + ret);
  ret
}

alpProp = function(a, g){ # device
  tmp = 0;
  Nalp = 0;

  for(n in 1:a$N)
    if(a$grp[n] != 2){
      tmp = tmp + a$y[n, g];
      Nalp = Nalp + 1;
    }
      
  avg = a$theAlp[a$m$theAlp] / a$gamAlp^2 + 
       tmp / (Nalp * a$sigAlp[a$m$sigAlp]^2)
  avg = avg / (1/a$gamAlp^2 + 1/a$sigAlp[a$m$sigAlp]^2)

  s = 1/sqrt(1/a$gamAlp^2 + 1/a$sigAlp[a$m$sigAlp]^2);

  u = runif(1);

#  if(u < piAlpPrime(a, g, avg, s)) {
  if(u < a$piAlp[a$m$piAlp]){
    new = 0;
  } else {
    new = sampleNormal(avg, s);
  }

  return(new);
}

sampleAlp_kernel1 = function(a){ # kernel <<<G, 1>>>
  for(g in 1:a$G){ 

    old = a$alp[a$m$alp, g];
    new = alpProp(a, g);
    
    lp = min(0, lAlp(a, g, new) - lAlp(a, g, old));
    lu = log(runif(1));
    
    if(lu < lp){ # accept
      a$alp[a$m$alp + 1, g] = new;
      a$acc$alp[g] = a$acc$alp[g] + 1;
    } else { # reject
      a$alp[a$m$alp + 1, g] = old;
    }
  }

  a
}

sampleAlp_kernel2 = function(a){ # kernel <<<1, 1>>>
  a$m$alp = a$m$alp + 1;
  a
}


sampleAlp = function(a){ # host
  a = sampleAlp_kernel1(a);
  a = sampleAlp_kernel2(a);
  
  return(a)
}

sampleTheAlp_kernel1 = function(a){ # kernel <<<G, 1>>>
  for(g in 1:a$G){ # PARALLELIZE
    if(a$alp[a$m$alp, g]){
      a$tmp1[g] = 1;
      a$tmp2[g] = a$alp[a$m$alp, g];
    } else {
      a$tmp1[g] = 0;
      a$tmp2[g] = 0;
    }
  }

  a
}

sampleTheAlp_kernel2 = function(a){ # parallel pairwise sum in Thrust
  Galp = 0;
  for(g in 1:a$G) # PARALLELIZE
    Galp = Galp + a$tmp1[g];

  a$s1 = Galp;
  a
}

sampleTheAlp_kernel3 = function(a){ # parallel pairwise sum in Thrust
  sm = 0;
  for(g in 1:a$G) # PARALLELIZE
    sm = sm + a$tmp2[g];

  a$s2 = sm;
  a
}

sampleTheAlp_kernel4 = function(a){ # kernel <<<1, 1>>>

  gs = a$gamAlp^2;
  ss = a$sigAlp[a$m$sigAlp]^2;
  den = (a$s1 * gs + ss);

  m = gs * a$s2 / den;
  s = gs * ss / den;

  a$theAlp[a$m$theAlp + 1] = sampleNormal(m, s);
  a$m$theAlp = a$m$theAlp + 1;

  a
}

sampleTheAlp = function(a){
  
  a = sampleTheAlp_kernel1(a);
  a = sampleTheAlp_kernel2(a);
  a = sampleTheAlp_kernel3(a);
  a = sampleTheAlp_kernel4(a);
 
  return(a);
}

sampleSigAlp_kernel1 = function(a){ # kernel <<<G, 1>>>

  for(g in 1:a$G){ # PARALLELIZE
    if(a$alp[a$m$alp, g]){
      a$tmp1[g] = (a$alp[a$m$alp, g] - a$theAlp[a$m$theAlp])^2;
      a$tmp2[g] = 1;
    } else {
      a$tmp1[g] = 0;
      a$tmp2[g] = 0;
    }
  }

  a
}

sampleSigAlp_kernel2 = function(a){ # parallel pairwise sum in Thrust
  rate = 0;  
  for(g in 1:a$G) # PARALLELIZE
    rate = rate + a$tmp1[g];

  a$s1 = rate;
  a
}

sampleSigAlp_kernel3 = function(a){ # parallel pairwise sum in Thrust
  Galp = 0;
  for(g in 1:a$G) # PARALLELIZE
    Galp = Galp + a$tmp2[g];  

  a$s2 = Galp;
  a
}


sampleSigAlp_kernel4 = function(a){ # parallel pairwise sum in Thrust
  shape = (a$s2 - 1) / 2;
  rate = a$s1 / 2;
  lb = 1/a$sigAlp0^2;

  if(shape >= 1 && rate > 0){
    a$sigAlp[a$m$sigAlp + 1] = 1/sqrt(sampleGamma(shape, rate, lb));
  } else {
    a$sigAlp[a$m$sigAlp + 1] = a$sigAlp[a$m$sigAlp]; 
  }

  a$m$sigAlp = a$m$sigAlp + 1;

  a;
}

sampleSigAlp = function(a){ # host

  a = sampleSigAlp_kernel1(a);
  a = sampleSigAlp_kernel2(a);
  a = sampleSigAlp_kernel3(a);
  a = sampleSigAlp_kernel4(a); 

  return(a)
}

samplePiAlp_kernel1 = function(a){ # kernel <<<1, 1>>>

  for(g in 1:a$G){ 
    if(a$alp[a$m$alp, g]){
      a$tmp1[g] = 1;
    } else {
      a$tmp1[g] = 0;
    }
  }

  a
}

samplePiAlp_kernel2 = function(a){ # pairwise sum in Thrust
  Galp = 0;
  for(g in 1:a$G) # PARALLELIZE
    Galp = Galp + a$tmp1[g];

  a$s1 = Galp;
  a  
}

samplePiAlp_kernel3 = function(a){ # kernel <<<1, 1>>>
  a$piAlp[a$m$piAlp + 1] = sampleBeta(a$G + a$s1 + a$aTau, a$s1 + a$bTau);
  a$m$piAlp = a$m$piAlp + 1;
  a
}

samplePiAlp = function(a){ # host

  a = samplePiAlp_kernel1(a);
  a = samplePiAlp_kernel2(a);  
  a = samplePiAlp_kernel3(a);

  return(a)
}

piDelPrime = function(a, g, avg, s){ # device
  prod = 1;

  for(n in 1:a$N){
    if(a$grp[n] == 2){
      tmp = a$c[a$m$c, n] + a$eps[a$m$eps, n, g] +
            mu(a, n, a$phi[a$m$phi, g], a$alp[a$m$del, g], 0);

      den = exp(a$y[n, g] * tmp - exp(tmp));

      i = a$c[a$m$c, n] + a$eps[a$m$eps, n, g] + a$phi[a$m$phi, g]

      A = -1/(2 * s * s) - 1/2;
      B = -i + avg/(s * s) + a$y[n, g] - 1;
      C = -(i * i)/2 + i * a$y[n, g] - i - (avg * avg)/(2 * s * s) - 1;
      D = 1.0/sqrt(2 * pi * s * s);
                     
      num = D * exp(C - (B * B)/(4 * A)) * sqrt(- pi/A);

      if(den != 0 && is.finite(num))
        prod = prod * num/den
    }
  }

  ret = ((1 - a$piAlp[a$m$piAlp]) / a$piAlp[a$m$piAlp]) * prod;
  ret = 1/(1 + ret);
  ret
}

delProp = function(a, g){ # device
  tmp = 0;
  Ndel = 0;

  for(n in 1:a$N)
    if(a$grp[n] != 2){
      tmp = tmp + a$y[n, g];
      Ndel = Ndel + 1;
    }
      
  avg = a$theDel[a$m$theDel] / a$gamDel^2 + 
       tmp / (Ndel * a$sigDel[a$m$sigDel]^2)
  avg = avg / (1/a$gamDel^2 + 1/a$sigDel[a$m$sigDel]^2)

  s = 1/sqrt(1/a$gamDel^2 + 1/a$sigDel[a$m$sigDel]^2);

  u = runif(1);

#  if(u < piDelPrime(a, g, avg, s)) {
  if(u < a$piDel[a$m$piDel]){
    new = 0;
  } else {
    new = sampleNormal(avg, s);
  }

  return(new)
}

sampleDel_kernel1 = function(a){ # kernel <<<G, 1>>>
  for(g in 1:a$G){ 

    old = a$del[a$m$del, g];
    new = delProp(a, g);
    
    lp = min(0, lDel(a, g, new) - lDel(a, g, old));
    lu = log(runif(1));
    
    if(lu < lp){ # accept
      a$del[a$m$del + 1, g] = new;
      a$acc$del[g] = a$acc$del[g] + 1;
    } else { # reject
      a$del[a$m$del + 1, g] = old;
    }
  }

  a
}

sampleDel_kernel2 = function(a){ # kernel <<<1 1>>>
  a$m$del = a$m$del + 1;
  a
}

sampleDel = function(a){ # host
  a = sampleDel_kernel1(a);
  a = sampleDel_kernel2(a);
  
  return(a)
}

sampleTheDel_kernel1 = function(a){ # kernel <<<G, 1>>>
  for(g in 1:a$G){ 
    if(a$del[a$m$del, g]){
      a$tmp1[g] = 1;
      a$tmp2[g] = a$del[a$m$del, g];
    } else {
      a$tmp1[g] = 0;
      a$tmp2[g] = 0;
    }
  }

  a
}

sampleTheDel_kernel2 = function(a){ # pairwise sum in Thrust
  Gdel = 0;
  for(g in 1:a$G) # PARALLELIZE  
    Gdel = Gdel + a$tmp1[g];

  a$s1 = Gdel;
  a
}

sampleTheDel_kernel3 = function(a){ # pairwise sum in Thrust
  sm = 0;
  for(g in 1:a$G) # PARALLELIZE
    sm = sm + a$tmp2[g];

  a$s2 = sm;
  a
}

sampleTheDel = function(a){ # host

  a = sampleTheDel_kernel1(a);
  a = sampleTheDel_kernel2(a);

  gs = a$gamDel^2;
  ss = a$sigDel[a$m$sigDel]^2;
  den = (a$s1 * gs + ss);

  m = gs * a$s2 / den;
  s = gs * ss / den;

  a$theDel[a$m$theDel + 1] = sampleNormal(m, s);
  a$m$theDel = a$m$theDel + 1;
  return(a)
}


sampleSigDel_kernel1 = function(a){ # kernel <<<G, 1>>>

  for(g in 1:a$G){ # PARALLELIZE
    if(a$del[a$m$del, g]){
      a$tmp1[g] = (a$del[a$m$del, g] - a$theDel[a$m$theDel])^2;
      a$tmp2[g] = 1;
    } else {
      a$tmp1[g] = 0;
      a$tmp2[g] = 0;
    }
  }
  
  a
}

sampleSigDel_kernel2 = function(a){ # pairwise sum in Thrust
  rate = 0;
  for(g in 1:a$G) # PARALLELIZE
    rate = rate + a$tmp1[g];

  a$s1 = rate;
  a
}

sampleSigDel_kernel3 = function(a){ # pairwise sum in Thrust
  Gdel = 0;
  for(g in 1:a$G) # PARALLELIZE
    Gdel = Gdel + a$tmp2[g];

  a$s2 = Gdel
  a
}

sampleSigDel_kernel4 = function(a){ # kernel <<<1, 1>>>
  shape = (a$s2 - 1) / 2;
  rate = a$s1 / 2;
  lb = 1/a$sigDel0^2;

  if(shape >= 1 && rate > 0){
    a$sigDel[a$m$sigDel + 1] = 1/sqrt(sampleGamma(shape, rate, lb));
  } else {
    a$sigDel[a$m$sigDel + 1] = a$sigDel[a$m$sigDel];
  }

  a$m$sigDel = a$m$sigDel + 1;

  a
}

sampleSigDel = function(a){ # host

  a = sampleSigDel_kernel1(a);
  a = sampleSigDel_kernel2(a);
  a = sampleSigDel_kernel3(a);
  a = sampleSigDel_kernel4(a);
  
  return(a)
}

samplePiDel_kernel1 = function(a){ # kernel <<<G, 1>>>

  for(g in 1:a$G){ # PARALLELIZE
    if(a$del[a$m$del, g]){
      a$tmp1[g] = 1; 
    } else {
      a$tmp1[g] = 0;
    }
  } 

  a
}

samplePiDel_kernel2 = function(a){ # pairwise sum in Thrust

  Gdel = 0;
  for(g in 1:a$G) # PARALLELIZE
     Gdel = Gdel + a$tmp1[g];

  a$s1 = Gdel;
  a
}

samplePiDel_kernel3 = function(a){ #kernel <<<1, 1>>>
  a$piDel[a$m$piDel + 1] = sampleBeta(a$G + a$s1 + a$aTau, a$s1 + a$bTau);
  a$m$piDel = a$m$piDel + 1;

  a
}

samplePiDel = function(a){ # host
  a = samplePiDel_kernel1(a);
  a = samplePiDel_kernel2(a);
  a = samplePiDel_kernel3(a);

  return(a)
}


samplePhiAlpDel_kernel1 = function(a){ # kernel <<<G, 1>>>
  for(g in 1:a$G){ # PARALLELIZE

    oldPhi = a$phi[a$m$phi, g];
    newPhi = sampleNormal(oldPhi, a$tun$phi[g]);

    oldAlp = a$alp[a$m$alp, g];
    newAlp = alpProp(a, g);

    oldDel = a$del[a$m$del, g];
    newDel = delProp(a, g);

    lp = min(0, lPhiAlpDel(a, g, newPhi, newAlp, newDel) 
              - lPhiAlpDel(a, g, oldPhi, oldAlp, oldDel));
    lu = log(runif(1));
    
    if(lu < lp){ # accept
      a$phi[a$m$phi + 1, g] = newPhi;
      a$alp[a$m$alp + 1, g] = newAlp;
      a$del[a$m$del + 1, g] = newDel;

      a$tun$phi[g] = a$tun$phi[g] * 1.1; 
      a$tun$alp[g] = a$tun$alp[g] * 1.1; 
      a$tun$del[g] = a$tun$del[g] * 1.1; 

      a$acc$phi[g] = a$acc$phi[g] + 1;
      a$acc$alp[g] = a$acc$alp[g] + 1;
      a$acc$del[g] = a$acc$del[g] + 1;
    } else { # reject
      a$phi[a$m$phi + 1, g] = oldPhi;
      a$alp[a$m$alp + 1, g] = oldAlp;
      a$del[a$m$del + 1, g] = oldDel;

      a$tun$phi[g] = a$tun$phi[g] / 1.1;
      a$tun$alp[g] = a$tun$alp[g] / 1.1;
      a$tun$del[g] = a$tun$del[g] / 1.1; 
    }
  }

  a
}

samplePhiAlpDel_kernel2 = function(a){ # kernel <<<1, 1>>>
  a$m$phi = a$m$phi + 1;
  a$m$alp = a$m$alp + 1;
  a$m$del = a$m$del + 1;
  a
}

samplePhiAlpDel = function(a){ # host
  a = samplePhiAlpDel_kernel1(a);
  a = samplePhiAlpDel_kernel2(a);
  return(a)
}

runChain = function(a, cfg){ # host
  for(m in 1:a$M){
    print(paste(m))

    a = sampleC(a);

    if(!cfg$skip$tau)
      a = sampleTau(a);

    if(!cfg$skip$piAlp)
      a = samplePiAlp(a);

    if(!cfg$skip$piDel)
      a = samplePiDel(a);

    if(!cfg$skip$d)
      a = sampleD(a);

    if(!cfg$skip$thePhi)
      a = sampleThePhi(a);

    if(!cfg$skip$theAlp)
      a = sampleTheAlp(a);

    if(!cfg$skip$theDel)
      a = sampleTheDel(a);

    if(!cfg$skip$sigC)
      a = sampleSigC(a); 

    if(!cfg$skip$sigPhi)
      a = sampleSigPhi(a); 

    if(!cfg$skip$sigAlp)
      a = sampleSigAlp(a);  

    if(!cfg$skip$sigDel)
      a = sampleSigDel(a);  

    a = sampleEta(a); 
    a = sampleEps(a);

    if(cfg$joint){
      a = samplePhiAlpDel(a);
    } else {
      a = samplePhi(a);
      a = sampleAlp(a);
      a = sampleDel(a);
    }
  }

  return(a)
}

allocSummary = function(a){ # device

  ret = list(
    sigC = rep(0, a$M),
    d = rep(0, a$M),
    tau = rep(0, a$M),

    thePhi = rep(0, a$M),
    sigPhi = rep(0, a$M),

    theAlp = rep(0, a$M),
    sigAlp = rep(0, a$M),
    piAlp = rep(0, a$M),

    theDel = rep(0, a$M),
    sigDel = rep(0, a$M),
    piDel = rep(0, a$M),

    acc = list(
      c = rep(0, a$N),
      eps = array(0, c(a$N, a$G)),
      d = 0,
      phi = rep(0, a$G),
      alp = rep(0, a$G),
      del = rep(0, a$G)
    )
  )

  return(ret)
}

summarizeChain = function(a){ # kernel <<<1, 1>>>
  ret = allocSummary(a);

  ret$sigC = a$sigC;
  ret$d = a$d;
  ret$tau = a$tau;

  ret$thePhi = a$thePhi;
  ret$sigPhi = a$sigPhi;

  ret$theAlp = a$theAlp;
  ret$sigAlp = a$sigAlp;
  ret$piAlp = a$piAlp;

  ret$theDel = a$theDel;
  ret$sigDel = a$sigDel;
  ret$piDel = a$piDel;

  ret$acc$d = a$acc$d / a$M;

  for(n in 1:a$N){
    ret$acc$c[n] = a$acc$c[n] / a$M;

    for(g in 1:a$G){
      ret$acc$eps[n, g] = a$acc$eps[n, g] / a$M;
    }
  }
  
  for(g in 1:a$G){
    ret$acc$phi[g] = a$acc$phi[g] / a$M;
    ret$acc$alp[g] = a$acc$alp[g] / a$M;
    ret$acc$del[g] = a$acc$del[g] / a$M;
  }

  ret
}

init = function(){
  h = hammer();
  y = h$y
  grp = h$grp
  newChain(y, grp, 10, 4, 100)
}


run = function(){
  cfg = getopts();
  a = newChain(cfg)
  a = runChain(a, cfg)
  s = summarizeChain(a)
  list(chain = a, summ = s)
}