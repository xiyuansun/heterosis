# Author: Will Landau
# PIs: Dan Nettleton, Jarad Niemi, Peng Liu
# Iowa State University
# May 2013
# 
# This code attempts to learn about differential expression
# between 2 treatment groups of an RNA-Seq dataset. See
# writeup.pdf for the model.

# sampling from known distributions

mySampleInt = function(N, n){

  # sample without replacement

  ret = rep(0, n);

  for(i in 1:n)
    ret[i] = -1;

  i = 1;
  while(i < n + 1){
    ret[i] = sample.int(N, 1);

    repeats = 0;
    j = 1;

    while(j < i){
      if(ret[j] == ret[i])
        repeats = repeats + 1;
 
      j = j + 1;
    }

    if(!repeats)
      i = i + 1;
  }

  # bubble sort

  if(n < 2)
    return(ret);

  swap = 1;
  while(swap){
    swap = 0;
    i = 2;

    while(i <= n){
      if(ret[i - 1] > ret[i]){
        tmp = ret[i]; 
        ret[i] = ret[i - 1];
        ret[i - 1] = tmp;
        
        swap = 1;
      }

      i = i + 1;
    }
  }

  ret
}

sampleNormal = function(m = 0, s = 1){ # host, device
  u1 = runif(1);
  u2 = runif(1);
  return(sqrt(-2 * log(u1)) * sin(2 * 3.14159265 * u2) * s + m);
}

sampleGamma = function(shape = 1, rate = 1, lb = 0){ # host, device
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

sampleBeta = function(a, b){ # host, device
  x = sampleGamma(a, 1, 0);
  y = sampleGamma(b, 1, 0);
  return(x / (x + y));
}

# data structures and initialization

allocConfig = function(){ # host
  list(
    datafile = "",
    groupfile = "",

    probsfile = "",
    hyperfile = "",
    ratesfile = "",
    parmsfile = "",

    probsFlag = 0,
    hyperFlag = 0,
    ratesFlag = 0,
    parmsFlag = 0,

    burnin = 0,
    joint = 0,
    heterosis = 1,

    M = 0,
    N = 0,
    G = 0,

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

    constSigC = 0,
    
    constD = 0,
    constTau = 0,
    
    constThePhi = 0,
    constTheAlp = 0,
    constTheDel = 0,

    constSigPhi = 0,
    constSigAlp = 0,
    constSigDel = 0,

    constPiAlp = 0,
    constPiDel = 0
      
  );
}

config = function(){ # host
  cfg = allocConfig();

  # default filenames and MCMC settings

  cfg$datafile = "../data/data.txt";
  cfg$groupfile = "../data/group.txt"; 

  cfg$probsfile = "../out/probs.txt";
  cfg$hyperfile = "../out/hyperparameters.txt";
  cfg$ratesfile = "../out/acceptance-rates.txt";
  cfg$parmsfile = "../out/example-parameters.txt";

  cfg$probsFlag = 1;
  cfg$hyperFlag = 1;
  cfg$ratesFlag = 1;
  cfg$parmsFlag = 1;

  cfg$burnin = 0;
  cfg$joint = 0;
  cfg$M = 10;

  # default values for initialization constants 

  cfg$sigC0 = 10; 
  cfg$d0 = 1e3; 

  cfg$aTau = 1e2; 
  cfg$aAlp = 1; 
  cfg$aDel = 1; 

  cfg$bTau = 1e2; 
  cfg$bAlp = 1; 
  cfg$bDel = 1; 

  cfg$gamPhi = 2; 
  cfg$gamAlp = 2; 
  cfg$gamDel = 2; 

  cfg$sigPhi0 = 2; 
  cfg$sigAlp0 = 2;
  cfg$sigDel0 = 2; 

  cfg$constSigC = 0;
  cfg$constD = 0;
  cfg$constTau = 0;

  cfg$constThePhi = 0;
  cfg$constTheAlp = 0;
  cfg$constTheDel = 0; 

  cfg$constSigPhi = 0; 
  cfg$constSigAlp = 0; 
  cfg$constSigDel = 0; 

  cfg$constPiAlp = 0;
  cfg$constPiDel = 0;

  # Use long getopt to parse command line args.

  # default initial values of hyperparameters

  if(!cfg$constSigC)
    cfg$sigC = runif(1, 0, cfg$sigC0);

  if(!cfg$constD)
    cfg$d = runif(1, 0, cfg$d0);

  if(!cfg$constTau)
    cfg$tau = sqrt(sampleGamma(shape = cfg$aTau, rate = cfg$bTau));

  if(!cfg$constThePhi)
    cfg$thePhi = sampleNormal(0, cfg$gamPhi);

  if(!cfg$constTheAlp)
    cfg$theAlp = sampleNormal(0, cfg$gamAlp);

  if(!cfg$constTheDel)
    cfg$theDel = sampleNormal(0, cfg$gamDel); 

  if(!cfg$constSigPhi)
    cfg$sigPhi = runif(1, 0, cfg$sigPhi0); 

  if(!cfg$constSigAlp)
    cfg$sigAlp = runif(1, 0, cfg$sigAlp0); 

  if(!cfg$constSigDel)
    cfg$sigDel = runif(1, 0, cfg$sigDel0); 

  if(!cfg$constPiAlp)
    cfg$piAlp = sampleBeta(cfg$aAlp, cfg$bAlp);

  if(!cfg$constPiDel)
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

    # current place in the chain for each parameter

    mC = 0,
    mSigC = 0,

    mEps = 0,
    mEta = 0,

    mD = 0,
    mTau = 0,

    mPhi = 0,
    mAlp = 0,
    mDel = 0,

    mThePhi = 0,
    mTheAlp = 0,
    mTheDel = 0,

    mSigPhi = 0,
    mSigAlp = 0,
    mSigDel = 0,

    mPiAlp = 0,
    mPiDel = 0,

    # tuning parameters for metropolis steps (std devs of normal dists)

    tuneC = rep(0, N),
    tuneEps = array(0, c(N, G)),
    tuneD = 0,
    tunePhi = rep(0, G),

    # number of acceptances for metropolis steps

    accD = 0,
    accC = rep(0, N),
    accEps = array(0, c(N, G)),
    accPhi = rep(0, G),
    accAlp = rep(0, G),
    accDel = rep(0, G)
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

  a$mC = 1;
  a$mSigC = 1;

  a$mEps = 1;
  a$mEta = 1;
  a$mD = 1;
  a$mTau = 1;

  a$mPhi = 1;
  a$mAlp = 1;
  a$mDel = 1;

  a$mThePhi = 1;
  a$mTheAlp = 1;
  a$mTheDel = 1;

  a$mSigPhi = 1;
  a$mSigAlp = 1;
  a$mSigDel = 1;

  a$mPiAlp = 1;
  a$mPiDel = 1;

  # tuning parameters for metropolis steps (std deviations of normal dists): 

  a$tuneD = 500;

  for(n in 1:a$N)
    a$tuneC[n] = 1;

  for(g in 1:a$G){
    a$tunePhi[g] = 1;

    for(n in 1:a$N)
      a$tuneEps[n, g] = 1;
  }

  # number of acceptances for metropolis steps

  a$accD = 0;

  for(n in 1:a$N){
    a$accC[n] = 0;
  
    for(g in 1:a$G)
      a$accEps[n, g] = 0;
  }

  for(g in 1:a$G){
    a$accPhi[g] = 0;
    a$accAlp[g] = 0;
    a$accDel[g] = 0;
  }
  
  a
}

readGrp = function(cfg, N){
  grp = scan(cfg$groupfile, quiet = T)

  unique = rep(0, N);
  nunique = 0;

  for(i in 1:N){
    
    match = 0;

    for(j in 1:max(1, nunique))
      if(grp[i] == unique[j])
        match = match + 1;

    if(!match){
      nunique = nunique + 1;
      unique[nunique] = grp[i];
    } 
  }

  if (nunique == 2){

    cfg$heterosis = 0;
    
    for(n in 1:N){
      if(grp[n] == unique[1]){
        grp[n] = 1;
      } else{
        grp[n] = 3;
      }
    }
  } else if (nunique != 3){
    print("ERROR: bad experimental design.")
  }

  list(grp = grp, cfg = cfg);
}

newChain = function(cfg){ # host (bunch of cudaMemcpies and kernels)

  y = t(as.matrix(read.table(cfg$datafile)))

  M = cfg$M;
  N = dim(y)[1];
  G = 50 # dim(y)[2];

  a = allocChain(M, N, G);

  l = readGrp(cfg, N);
  cfg = l$cfg;
  grp = l$grp;

  cfg$N = N;
  cfg$G = G;

  a$heterosis = cfg$heterosis; # CudaMemcpy
  a$parmsFlag = cfg$parmsFlag; # CudaMemcpy

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

  list(a = a, cfg = cfg);
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
    a$tmp1[g] = exp(a$eps[a$mEps, n, g] + mu(a, n, a$phi[a$mPhi, g], 
                                            a$alp[a$mAlp, g], a$del[a$mDel, g]));

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
        (2 * a$sigC[a$mSigC] * a$sigC[a$mSigC]);

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
        exp(a$c[a$mC, n] + arg + mu(a, n, a$phi[a$mPhi, g], 
                                     a$alp[a$mAlp, g], a$del[a$mDel, g])) 
          - arg^2 / (2 * a$eta[a$mEta, g]^2);

  return(ret);
}

lD_kernel1 = function(a){ # kernel <<<G, 1>>>
  for(g in 1:a$G){ 
    a$tmp1[g] = 2 * log(a$eta[a$mEta, g]);
    a$tmp2[g] = 1/(a$eta[a$mEta, g] * a$eta[a$mEta, g]);
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

  a$tmp1[1] = arg * a$tau[a$mTau]^2 / 2;

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
    tmp = mu(a, n, arg, a$alp[a$mAlp, g], a$del[a$mDel, g]);
    s = s + a$y[n, g] * tmp - exp(a$c[a$mC, n] + 
        a$eps[a$mEps, n, g] + tmp);
  }
 
  ret = s - (arg - a$thePhi[a$mThePhi])^2 / (2 * a$sigPhi[a$mSigPhi]^2);
  return(ret);
}

lAlp = function(a, g, arg){ # device
  
  s = 0; 
  for(n in 1:a$N){
    if(a$grp[n] != 2){
      tmp = mu(a, n, a$phi[a$mPhi, g], arg, a$del[a$mDel, g]);
      s = s + a$y[n, g] * tmp - exp(a$c[a$mC, n] + 
          a$eps[a$mEps, n, g] + tmp);
    }
  }
 
  if(arg * arg > 1e-6){
    tmp = -(arg - a$theAlp[a$mTheAlp])^2 / (2 * a$sigAlp[a$mSigAlp]^2) -
                log(1 - a$piAlp[a$mPiAlp]);
  } else {
    tmp = log(a$piAlp[a$mPiAlp]);
  }

  ret = s + tmp;
  return(ret);
}

lDel = function(a, g, arg){ # device 
  
  s = 0; 
  for(n in 1:a$N){
    if(a$grp[n] != 2){
      tmp = mu(a, n, a$phi[a$mPhi, g], a$alp[a$mAlp, g], arg);
      s = s + a$y[n, g] * tmp - exp(a$c[a$mC, n] + 
          a$eps[a$mEps, n, g] + tmp);
    }
  }
 
  if(arg * arg > 1e-6){
    tmp = -(arg - a$theDel[a$mTheDel])^2 / (2 * a$sigDel[a$mSigDel]^2) -
                log(1 - a$piDel[a$mPiDel]);
  } else {
    tmp = log(a$piDel[a$mPiDel]);
  }

  ret = s + tmp;
  return(ret);
}

lPhiAlpDel = function(a, g, argPhi, argAlp, argDel){ # device
 
  s = 0; 
  tmp = 0;

  for(n in 1:a$N){
    tmp = mu(a, n, argPhi, argAlp, argDel);
    s = s + a$y[n, g] * tmp - exp(a$c[a$mC, n] + 
        a$eps[a$mEps, n, g] + tmp);
  }

  # phi part 
  ret = s - (argPhi - a$thePhi[a$mThePhi])^2 / (2 * a$sigPhi[a$mSigPhi]^2);

  # alpha part
  if(argAlp * argAlp > 1e-6){
    tmp = -(argAlp - a$theAlp[a$mTheAlp])^2 / (2 * a$sigAlp[a$mSigAlp]^2) -
                log(1 - a$piAlp[a$mPiAlp]);
  } else {
    tmp = log(a$piAlp[a$mPiAlp]);
  }

  ret = ret + tmp;

  # delta part
  if(argDel * argDel > 1e-6){
    tmp = -(argDel - a$theDel[a$mTheDel])^2 / (2 * a$sigDel[a$mSigDel]^2) -
                log(1 - a$piDel[a$mPiDel]);
  } else {
    tmp = log(a$piDel[a$mPiDel]);
  }

  ret = ret + tmp

  return(ret);
}



# samplers

sampleC_kernel1 = function(a, n){ # kernel <<<1, 1>>>
  a$old[1, 1] = a$c[a$mC, n];
  a$new[1, 1] = sampleNormal(a$old[1, 1], a$tuneC[n]);

  a
}

sampleC_kernel2 = function(a, n){ # kernel <<<1, 1>>>
  lp = min(0, a$lNew[1, 1] - a$lOld[1, 1]);
  lu = log(runif(1));
    
  if(lu < lp){ # accept
    a$c[a$mC + 1, n] = a$new[1, 1];
    a$tuneC[n] = a$tuneC[n] * 1.1; # Increase the proposal variance to avoid  
                                   # getting stuck in a mode
    a$accC[n] = a$accC[n] + 1;
  } else { # reject
    a$c[a$mC + 1, n] = a$old[1, 1];
    a$tuneC[n] = a$tuneC[n] / 1.1; # If you're rejecting too often, decrease the  
                                   # proposal variance to sample closer to 
                                   # the last accepted value.
  }

  a
}

sampleC_kernel3 = function(a){ # kernel <<<1, 1>>>
  a$mC = a$mC + 1;
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
    a$tmp1[n] = a$c[a$mC, n]^2;

  rate = 0;
  for(n in 1:a$N) 
    rate = rate + a$tmp1[n];
  
  shape = (a$N - 1) / 2; 
  rate = rate / 2;
  lb = 1 / a$sigC0^2  

  if(shape >= 1 && rate > 0){
    a$sigC[a$mSigC + 1] = 1/sqrt(sampleGamma(shape, rate, lb));
  } else {
    a$sigC[a$mSigC + 1] = a$sigC[a$mSigC];
  }

  a$mSigC = a$mSigC + 1;
  return(a)
}

sampleEps_kernel1 = function(a){ # kernel <<<N, G>>>
  for(g in 1:a$G){
    for(n in 1:a$N){ 
      old = a$eps[a$mEps, n, g];
      new = sampleNormal(old, a$tuneEps[n, g]);

      lp = min(0, lEps(a, n, g, new) - lEps(a, n, g, old));
      lu = log(runif(1));
      
      if(lu < lp){ # accept
        a$eps[a$mEps + 1, n, g] = new;
        a$tuneEps[n, g] = a$tuneEps[n, g] * 1.1;
        a$accEps[n, g] = a$accEps[n, g] + 1; 
      } else { # reject
        a$eps[a$mEps + 1, n, g] = old;
        a$tuneEps[n, g] = a$tuneEps[n, g] / 1.1;
      }
    }
  }

  a
}

sampleEps_kernel2 = function(a){ # kernel <<<1, 1>>>
  a$mEps = a$mEps + 1;
  a
}

sampleEps = function(a){ # host

  a = sampleEps_kernel1(a);
  a = sampleEps_kernel2(a);

  a
}

sampleEta_kernel1 = function(a){ # kernel <<<1, 1>>>
  a$shape = (a$N + a$d[a$mD]) / 2; 
  a
}

sampleEta_kernel2 = function(a){ # kernel <<<G, 1>>>
  for(g in 1:a$G){

    for(n in 1:a$N) 
      a$tmp1[n] = a$eps[a$mEps, n, g]^2;

    rate = 0;
    for(n in 1:a$N) 
      rate = rate + a$tmp1[n];
  
    rate = (rate + a$d[a$mD] * a$tau[a$mTau] * a$tau[a$mTau]) / 2; 

    if(a$shape >= 1 && rate > 0){
      a$eta[a$mEta + 1, g] = 1/sqrt(sampleGamma(a$shape, rate));
    } else {
      a$eta[a$mEta + 1, g] = a$eta[a$mEta, g];
    }
  }
  a
}

sampleEta_kernel3 = function(a){ # kernel <<<1, 1>>>
  a$mEta = a$mEta + 1;
  a
}

sampleEta = function(a){

  a = sampleEta_kernel1(a);
  a = sampleEta_kernel2(a);
  a = sampleEta_kernel3(a);

  return(a)
}


sampleD_kernel1 = function(a){ # kernel <<<1, 1>>>
  a$old[1, 1] = a$d[a$mD];
  a$new[1, 1] = sampleNormal(a$old[1, 1], a$tuneD);
  a
}

sampleD_kernel2 = function(a){ # kernel <<<1, 1>>>
  lp = min(0, a$lNew[1, 1] - a$lOld[1, 1]);
  lu = log(runif(1));

  if(lu < lp){ # accept
    a$d[a$mD + 1] = a$new[1, 1];
    a$tuneD = a$tuneD * 1.1; # Increase the proposal variance to avoid getting 
                                 # stuck in a mode
    a$accD = a$accD + 1;
  } else { # reject
    a$d[a$mD + 1] = a$old[1, 1];
    a$tuneD = a$tuneD / 1.1; # If you're rejecting too often, decrease the proposal 
                                 # variance to sample closer to the last accepted value.
  }

  a$mD = a$mD + 1;
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
    a$tmp1[g] = 1/a$eta[a$mEta, g]^2;

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
  rate = a$rate * a$d[a$mD] / 2 + a$bTau;
  shape = a$aTau + a$G * a$d[a$mD] / 2;

  if(shape >= 1 && rate > 0){
    a$tau[a$mTau + 1] = 1/sqrt(sampleGamma(shape, rate));
  } else {
    a$tau[a$mTau + 1] = a$tau[a$mTau];
  }

  a$mTau = a$mTau + 1;

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

    old = a$phi[a$mPhi, g];
    new = sampleNormal(old, a$tunePhi[g]);

    lp = min(0, lPhi(a, g, new) - lPhi(a, g, old));
    lu = log(runif(1));
    
    if(lu < lp){ # accept
      a$phi[a$mPhi + 1, g] = new;
      a$tunePhi[g] = a$tunePhi[g] * 1.1; 
      a$accPhi[g] = a$accPhi[g] + 1;
    } else { # reject
      a$phi[a$mPhi + 1, g] = old;
      a$tunePhi[g] = a$tunePhi[g] / 1.1; 
    }
  }

  a
}

samplePhi_kernel2 = function(a){ # kernel <<<1, 1>>>
  a$mPhi = a$mPhi + 1;
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
    a$s1 = a$s1 + a$phi[a$mPhi, g];

  a
}

sampleThePhi_kernel2 = function(a){ # kernel <<<1, 1>>>
  gs = a$gamPhi^2;
  ss = a$sigPhi[a$mSigPhi]^2;
  den = (a$G * gs + ss);

  m = gs * a$s1 / den;
  s = gs * ss / den;

  a$thePhi[a$mThePhi + 1] = sampleNormal(m, s);
  a$mThePhi = a$mThePhi + 1;

  a
}

sampleThePhi = function(a){ # host
  a = sampleThePhi_kernel1(a);
  a = sampleThePhi_kernel2(a);
  return(a)
}

sampleSigPhi_kernel1 = function(a){ # kernel <<<G, 1>>>
  for(g in 1:a$G) # PARALLELIZE
    a$tmp1[g] = (a$phi[a$mPhi, g] - a$thePhi[a$mThePhi])^2;

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
    a$sigPhi[a$mSigPhi + 1] = 1/sqrt(sampleGamma(shape, rate, lb));
  } else {
    a$sigPhi[a$mSigPhi + 1] = a$sigPhi[a$mSigPhi];
  }

  a$mSigPhi = a$mSigPhi + 1;

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
      tmp = a$c[a$mC, n] + a$eps[a$mEps, n, g] +
            mu(a, n, a$phi[a$mPhi, g], 0, a$del[a$mDel, g]);

      den = exp(a$y[n, g] * tmp - exp(tmp));

      i = a$c[a$mC, n] + a$eps[a$mEps, n, g] + a$phi[a$mPhi, g]

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

  ret = ((1 - a$piAlp[a$mPiAlp]) / a$piAlp[a$mPiAlp]) * prod;
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
      
  avg = a$theAlp[a$mTheAlp] / a$gamAlp^2 + 
       tmp / (Nalp * a$sigAlp[a$mSigAlp]^2)
  avg = avg / (1/a$gamAlp^2 + 1/a$sigAlp[a$mSigAlp]^2)

  s = 1/sqrt(1/a$gamAlp^2 + 1/a$sigAlp[a$mSigAlp]^2);

  u = runif(1);

#  if(u < piAlpPrime(a, g, avg, s)) {
  if(u < a$piAlp[a$mPiAlp]){
    new = 0;
  } else {
    new = sampleNormal(avg, s);
  }

  return(new);
}

sampleAlp_kernel1 = function(a){ # kernel <<<G, 1>>>
  for(g in 1:a$G){ 

    old = a$alp[a$mAlp, g];
    new = alpProp(a, g);
    
    lp = min(0, lAlp(a, g, new) - lAlp(a, g, old));
    lu = log(runif(1));
    
    if(lu < lp){ # accept
      a$alp[a$mAlp + 1, g] = new;
      a$accAlp[g] = a$accAlp[g] + 1;
    } else { # reject
      a$alp[a$mAlp + 1, g] = old;
    }
  }

  a
}

sampleAlp_kernel2 = function(a){ # kernel <<<1, 1>>>
  a$mAlp = a$mAlp + 1;
  a
}


sampleAlp = function(a){ # host
  a = sampleAlp_kernel1(a);
  a = sampleAlp_kernel2(a);
  
  return(a)
}

sampleTheAlp_kernel1 = function(a){ # kernel <<<G, 1>>>
  for(g in 1:a$G){ # PARALLELIZE
    if(a$alp[a$mAlp, g]){
      a$tmp1[g] = 1;
      a$tmp2[g] = a$alp[a$mAlp, g];
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
  ss = a$sigAlp[a$mSigAlp]^2;
  den = (a$s1 * gs + ss);

  m = gs * a$s2 / den;
  s = gs * ss / den;

  a$theAlp[a$mTheAlp + 1] = sampleNormal(m, s);
  a$mTheAlp = a$mTheAlp + 1;

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
    if(a$alp[a$mAlp, g]){
      a$tmp1[g] = (a$alp[a$mAlp, g] - a$theAlp[a$mTheAlp])^2;
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
    a$sigAlp[a$mSigAlp + 1] = 1/sqrt(sampleGamma(shape, rate, lb));
  } else {
    a$sigAlp[a$mSigAlp + 1] = a$sigAlp[a$mSigAlp]; 
  }

  a$mSigAlp = a$mSigAlp + 1;

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
    if(a$alp[a$mAlp, g]){
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
  a$piAlp[a$mPiAlp + 1] = sampleBeta(a$G + a$s1 + a$aTau, a$s1 + a$bTau);
  a$mPiAlp = a$mPiAlp + 1;
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
      tmp = a$c[a$mC, n] + a$eps[a$mEps, n, g] +
            mu(a, n, a$phi[a$mPhi, g], a$alp[a$mDel, g], 0);

      den = exp(a$y[n, g] * tmp - exp(tmp));

      i = a$c[a$mC, n] + a$eps[a$mEps, n, g] + a$phi[a$mPhi, g]

      A = -1/(2 * s * s) - 1/2;
      B = -i + avg/(s * s) + a$y[n, g] - 1;
      C = -(i * i)/2 + i * a$y[n, g] - i - (avg * avg)/(2 * s * s) - 1;
      D = 1.0/sqrt(2 * pi * s * s);
                     
      num = D * exp(C - (B * B)/(4 * A)) * sqrt(- pi/A);

      if(den != 0 && is.finite(num))
        prod = prod * num/den
    }
  }

  ret = ((1 - a$piAlp[a$mPiAlp]) / a$piAlp[a$mPiAlp]) * prod;
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
      
  avg = a$theDel[a$mTheDel] / a$gamDel^2 + 
       tmp / (Ndel * a$sigDel[a$mSigDel]^2)
  avg = avg / (1/a$gamDel^2 + 1/a$sigDel[a$mSigDel]^2)

  s = 1/sqrt(1/a$gamDel^2 + 1/a$sigDel[a$mSigDel]^2);

  u = runif(1);

#  if(u < piDelPrime(a, g, avg, s)) {
  if(u < a$piDel[a$mPiDel]){
    new = 0;
  } else {
    new = sampleNormal(avg, s);
  }

  return(new)
}

sampleDel_kernel1 = function(a){ # kernel <<<G, 1>>>
  for(g in 1:a$G){ 

    old = a$del[a$mDel, g];
    new = delProp(a, g);
    
    lp = min(0, lDel(a, g, new) - lDel(a, g, old));
    lu = log(runif(1));
    
    if(lu < lp){ # accept
      a$del[a$mDel + 1, g] = new;
      a$accDel[g] = a$accDel[g] + 1;
    } else { # reject
      a$del[a$mDel + 1, g] = old;
    }
  }

  a
}

sampleDel_kernel2 = function(a){ # kernel <<<1 1>>>
  a$mDel = a$mDel + 1;
  a
}

sampleDel = function(a){ # host
  a = sampleDel_kernel1(a);
  a = sampleDel_kernel2(a);
  
  return(a)
}

sampleTheDel_kernel1 = function(a){ # kernel <<<G, 1>>>
  for(g in 1:a$G){ 
    if(a$del[a$mDel, g]){
      a$tmp1[g] = 1;
      a$tmp2[g] = a$del[a$mDel, g];
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
  ss = a$sigDel[a$mSigDel]^2;
  den = (a$s1 * gs + ss);

  m = gs * a$s2 / den;
  s = gs * ss / den;

  a$theDel[a$mTheDel + 1] = sampleNormal(m, s);
  a$mTheDel = a$mTheDel + 1;
  return(a)
}


sampleSigDel_kernel1 = function(a){ # kernel <<<G, 1>>>

  for(g in 1:a$G){ # PARALLELIZE
    if(a$del[a$mDel, g]){
      a$tmp1[g] = (a$del[a$mDel, g] - a$theDel[a$mTheDel])^2;
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
    a$sigDel[a$mSigDel + 1] = 1/sqrt(sampleGamma(shape, rate, lb));
  } else {
    a$sigDel[a$mSigDel + 1] = a$sigDel[a$mSigDel];
  }

  a$mSigDel = a$mSigDel + 1;

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
    if(a$del[a$mDel, g]){
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
  a$piDel[a$mPiDel + 1] = sampleBeta(a$G + a$s1 + a$aTau, a$s1 + a$bTau);
  a$mPiDel = a$mPiDel + 1;

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

    oldPhi = a$phi[a$mPhi, g];
    newPhi = sampleNormal(oldPhi, a$tunePhi[g]);

    oldAlp = a$alp[a$mAlp, g];
    newAlp = alpProp(a, g);

    oldDel = a$del[a$mDel, g];
    newDel = delProp(a, g);

    lp = min(0, lPhiAlpDel(a, g, newPhi, newAlp, newDel) 
              - lPhiAlpDel(a, g, oldPhi, oldAlp, oldDel));
    lu = log(runif(1));
    
    if(lu < lp){ # accept
      a$phi[a$mPhi + 1, g] = newPhi;
      a$alp[a$mAlp + 1, g] = newAlp;
      a$del[a$mDel + 1, g] = newDel;

      a$tunePhi[g] = a$tunePhi[g] * 1.1; 

      a$accPhi[g] = a$accPhi[g] + 1;
      a$accAlp[g] = a$accAlp[g] + 1;
      a$accDel[g] = a$accDel[g] + 1;
    } else { # reject
      a$phi[a$mPhi + 1, g] = oldPhi;
      a$alp[a$mAlp + 1, g] = oldAlp;
      a$del[a$mDel + 1, g] = oldDel;

      a$tunePhi[g] = a$tunePhi[g] / 1.1;
    }
  }

  a
}

samplePhiAlpDel_kernel2 = function(a){ # kernel <<<1, 1>>>
  a$mPhi = a$mPhi + 1;
  a$mAlp = a$mAlp + 1;
  a$mDel = a$mDel + 1;
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

    if(!cfg$constTau)
      a = sampleTau(a);

    if(!cfg$constPiAlp)
      a = samplePiAlp(a);

    if(!cfg$constPiDel && cfg$heterosis)
      a = samplePiDel(a);

    if(!cfg$constD)
      a = sampleD(a);

    if(!cfg$constThePhi)
      a = sampleThePhi(a);

    if(!cfg$constTheAlp)
      a = sampleTheAlp(a);

    if(!cfg$constTheDel && cfg$heterosis)
      a = sampleTheDel(a);

    if(!cfg$constSigC)
      a = sampleSigC(a); 

    if(!cfg$constSigPhi)
      a = sampleSigPhi(a); 

    if(!cfg$constSigAlp)
      a = sampleSigAlp(a);  

    if(!cfg$constSigDel && cfg$heterosis)
      a = sampleSigDel(a);  

    a = sampleEta(a); 
    a = sampleEps(a);

    if(cfg$joint && cfg$heterosis){
      a = samplePhiAlpDel(a);
    } else {
      a = samplePhi(a);
      a = sampleAlp(a);

      if(cfg$heterosis)
        a = sampleDel(a);
    }
  }

  return(a)
}

allocSummary = function(a){ # host, device

  ret = list(

    # hyperparamters

    sigC = rep(0, a$M + 1),
    d = rep(0, a$M + 1),
    tau = rep(0, a$M + 1),

    thePhi = rep(0, a$M + 1),
    theAlp = rep(0, a$M + 1),
    theDel = rep(0, a$M + 1),

    sigPhi = rep(0, a$M + 1),
    sigAlp = rep(0, a$M + 1),
    sigDel = rep(0, a$M + 1),

    piAlp = rep(0, a$M + 1),
    piDel = rep(0, a$M + 1),

    # acceptance rates of metropolis steps

    accD = 0,
    accC = rep(0, a$N),
    accEps = rep(0, a$G),
    accPhi = rep(0, a$G),
    accAlp = rep(0, a$G),
    accDel = rep(0, a$G),

    # probabilities of differential expression and heterosis

    prob_de = rep(0, a$G),
    prob_hph = rep(0, a$G),
    prob_lph = rep(0, a$G),
    prob_mph = rep(0, a$G),

    # samples of parameters from 5 random example libraries 
    #   and 5 random example genes

    libs = rep(0, 5),
    genes = rep(0, 5),

    c = array(0, c(a$M + 1, 5)),
    eps = array(0, c(a$M + 1, 5, 5)),
    eta = array(0, c(a$M + 1, 5)),
    phi = array(0, c(a$M + 1, 5)),
    alp = array(0, c(a$M + 1, 5)),
    del = array(0, c(a$M + 1, 5))
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

  ret$accD = a$accD / a$M;

  for(n in 1:a$N)
    ret$accC[n] = a$accC[n] / a$M;
  
  for(g in 1:a$G){

    for(n in 1:a$N)
      ret$accEps[g] = ret$accEps[g] + a$accEps[n, g];

    ret$accEps[g] = ret$accEps[g] / (a$M * a$N);

    ret$accPhi[g] = a$accPhi[g] / a$M;
    ret$accAlp[g] = a$accAlp[g] / a$M;
    ret$accDel[g] = a$accDel[g] / a$M;

    ret$prob_de[g] = 0;
    ret$prob_hph[g] = 0;
    ret$prob_lph[g] = 0;
    ret$prob_mph[g] = 0;
  }

  for(m in (2 + a$burnin):(a$M + 1)){
    for(g in 1:a$G){

      if(a$alp[m, g] > 1e-6)
        ret$prob_de[g] = ret$prob_de[g] + 1;

      if(a$heterosis){

        if(a$del[m, g] > sqrt(a$alp[m, g] * a$alp[m, g]))
          ret$prob_hph[g] = ret$prob_hph[g] + 1;

        if(a$del[m, g] < -sqrt(a$alp[m, g] * a$alp[m, g]))
          ret$prob_lph[g] = ret$prob_lph[g] + 1;

        if(a$del[m, g] > 1e-6)
          ret$prob_mph[g] = ret$prob_mph[g] + 1;
      }
    }
  }

  for(g in 1:a$G){
    ret$prob_de[g] = ret$prob_de[g] / (a$M - a$burnin);

    if(a$heterosis){
      ret$prob_hph[g] = ret$prob_hph[g] / (a$M - a$burnin);
      ret$prob_lph[g] = ret$prob_lph[g] / (a$M - a$burnin);
      ret$prob_mph[g] = ret$prob_mph[g] / (a$M - a$burnin);
    }
  }

  if(a$parmsFlag){

    libs = mySampleInt(a$N, 5);
    genes = mySampleInt(a$G, 5);

    for(i in 1:5){
      ret$libs[i] = libs[i];
      ret$genes[i] = genes[i];
    }

    for(m in 1:(a$M + 1)){
      for(n in 1:5){
        ret$c[m, n] = a$c[m, libs[n]];

        for(g in 1:5){
          ret$eps[m, n, g] = a$eps[m, libs[n], genes[g]];
        }
      }

      for(g in 1:5){
        ret$eta[m, g] = a$eta[m, genes[g]]; 
        ret$phi[m, g] = a$phi[m, genes[g]];
        ret$alp[m, g] = a$alp[m, genes[g]];
        ret$del[m, g] = a$del[m, genes[g]];
      }
    }
  }

  ret
}

printSummary = function(s, cfg){

  if(cfg$probsFlag){

    str = "de"
    if(cfg$heterosis)
      str = paste(str, "hph", "lph", "mph")

    if(cfg$heterosis){
      for(g in 1:cfg$G)
        str = c(str, paste(c(s$prob_de[g], s$prob_hph[g], s$prob_lph[g], 
                         s$prob_mph[g]), collapse = " "))
    } else {
      for(g in 1:cfg$G)
        str = c(str, paste(s$prob_de[g]))
    }

    fp = file(cfg$probsfile);
    writeLines(str, fp)
    close(fp);
  }

  if(cfg$hyperFlag){

    str = "sigC d tau thePhi theAlp theDel sigPhi sigAlp sigDel piAlp piDel"

    for(m in 1:(cfg$M + 1))
      str = c(str, paste(round(c(s$sigC[m], s$d[m], s$tau[m], s$thePhi[m], 
                         s$theAlp[m], s$theDel[m], s$sigPhi[m], s$sigAlp[m],  
                         s$sigDel[m], s$piAlp[m], s$piDel[m]), 2), 
              collapse = " "))

    fp = file(cfg$hyperfile);
    writeLines(str, fp)
    close(fp);

  }

  if(cfg$ratesFlag){

    str = "d c meanEps phi alp del";

    for(g in 1:cfg$G){
      line = c();
      
      if(g == 1){
        line = c(line, s$accD);
      } else {
        line = c(line, 0/0);
      }

      if(g <= cfg$N){
        line = c(line, s$accC[g]);
      } else {
        line = c(line, 0/0);
      }

      line = c(line, s$accEps[g], s$accPhi[g], s$accAlp[g], s$accDel[g]);
      line = round(line, 2);
      line = paste(line, collapse = " ");

      str = c(str, line);
    }

    fp = file(cfg$ratesfile);
    writeLines(str, fp)
    close(fp);

  }

  if(cfg$parmsFlag){
    
    str = paste("c", s$libs, sep = "");
    str = c(str, paste(c("phi"), s$genes, sep = ""));
    str = c(str, paste(c("alp"), s$genes, sep = ""));
    str = c(str, paste(c("del"), s$genes, sep = ""));
    str = c(str, paste(c("eta"), s$genes, sep = ""));  

    for(n in 1:5)
      for(g in 1:5)
        str = c(str, paste("eps_lib", s$libs[n], "_gene", s$genes[g], sep=""));
   
    str = paste(str, collapse = " ");

    for(m in 1:(cfg$M + 1)){
      line = c();

      for(n in 1:5)
        line = c(line, s$c[m, n]);

      for(g in 1:5)
        line = c(line, s$phi[m, g]);

      for(g in 1:5)
        line = c(line, s$alp[m, g]);

      for(g in 1:5)
        line = c(line, s$del[m, g]);

      for(g in 1:5)
        line = c(line, s$eta[m, g]);

      for(n in 1:5)
        for(g in 1:5)
          line = c(line, s$eps[m, n, g]);
 
      line = round(line, 2);
      line = paste(line, collapse = " ");
      str = c(str, line);
    }
    
    fp = file(cfg$parmsfile);
    writeLines(str, fp)
    close(fp);
  }
}


run = function(){
  cfg = config();

  l = newChain(cfg);
  cfg = l$cfg;
  a = l$a;

  a = runChain(a, cfg);
  s = summarizeChain(a);
  printSummary(s, cfg);
}