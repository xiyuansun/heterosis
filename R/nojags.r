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

# reading data

hammer = function(){ # host
  load(url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/hammer_eset.RData"))
  counts = exprs(hammer.eset)
  counts = counts[rowSums(counts) > 0,]
  
  group = as.vector(phenoData(hammer.eset)$protocol)
  group[group == "control"] = "1"
  group[group == "L5 SNL"] = "2"
  group = as.factor(as.numeric(group))

  return(list(y = t(counts), 
       grp = group, 
       G = dim(counts)[1], 
       N = dim(counts)[2]))
}

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

  if(shape >= 1){ # Marsaglia and Tsang (2000)

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

allocChain = function(y, M, N, G){ # host (bunch of cudaMallocs)
 list(

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

    # data 

    y = array(0, c(N, G)),
    yMeanG = rep(0, N),
    grp = rep(0, N),

    M = 0,
    N = 0,
    G = 0,

    # parameters

    c = array(0, c(M, N)),
      sigC = rep(0, M),
    
    eps = array(0, c(M, N, G)),
      eta = array(0, c(M, G)),
        d = rep(0, M),
        tau = rep(0, M),

    phi = array(0, c(M, G)),
      thePhi = rep(0, M),
      sigPhi = rep(0, M),

    alp = array(0, c(M, G)),
      theAlp = rep(0, M),
      sigAlp = rep(0, M),
      piAlp = rep(0, M),

    del = array(0, c(M, G)),
      theDel = rep(0, M),
      sigDel = rep(0, M),
      piDel = rep(0, M),

    # temporary and return values

    tmp1 = rep(0, G),
    tmp2 = rep(0, G),  

    shape = 0,
    rate = 0,
    
    old = 0,
    new = 0,

    lOld = 0,
    lNew = 0,

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
    )
  );
}

newChain_kernel1 = function(a){ # kernel: 1 block, 1 thread each
  a$sigC[1] = runif(1, 0, a$sigC0)

  a$d[1] = runif(1, 0, a$d0)
  a$tau[1] = sqrt(sampleGamma(shape = a$aTau, rate = a$bTau))

  a$thePhi[1] = sampleNormal(0, a$gamPhi)
  a$theAlp[1] = sampleNormal(0, a$gamAlp)
  a$theDel[1] = sampleNormal(0, a$gamDel)

  a$sigPhi[1] = runif(1, 0, a$sigPhi0)
  a$sigAlp[1] = runif(1, 0, a$sigAlp0)
  a$sigDel[1] = runif(1, 0, a$sigDel0)

  a$piAlp[1] = sampleBeta(a$aAlp, a$bAlp)
  a$piDel[1] = sampleBeta(a$aDel, a$bDel)

  return(a);
}

newChain_kernel2 = function(a){ # kernel: G blocks, 1 thread each
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

newChain = function(y, grp, M, N, G){ # host (bunch of cudaMemCpies and kernels)

  a = allocChain(y, M, N, G);

  # data: cudaMemCpies

  for(n in 1:N){
    a$grp[n] = grp[n];
    tmp = 0;

    for(g in 1:G){
      a$y[n, g] = y[n, g];
      tmp = tmp + y[n, g];
    }

    a$yMeanG[n] = tmp / G;
  }

  a$M = M;
  a$N = N;
  a$G = G;

  # initialization constants: cudaMemCpies

  a$sigC0 = 10;
  a$d0 = 1e3;

  a$aTau = 1e2;
  a$aAlp = 1;
  a$aDel = 1;

  a$bTau = 1e2;
  a$bAlp = 1;
  a$bDel = 1;
  
  a$gamPhi = 2;
  a$gamAlp = 2;
  a$gamDel = 2;

  a$sigPhi0 = 2;
  a$sigAlp0 = 2;
  a$sigDel0 = 2;

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
    a$c[1, n] = lqts[n] - s; # cudaMemCpy

  # location in chain: cudaMemCpies

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

  # tuning parameters for metropolis steps (std deviations of normal dists): cudaMemCpies

  for(n in 1:N)
    a$tun$c[n] = 1;

  for(g in 1:G){
    a$tun$phi[g] = 1;

    for(n in 1:N)
      a$tun$eps[n, g] = 1;
  }

  a$tun$d = 500;

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
    arg = a$new;
  } else {
    arg = a$old;
  }

  ret = arg * a$G * a$yMeanG[n] - exp(arg) * a$tmp2[1] - (arg*arg) / 
        (2 * a$sigC[a$m$sigC] * a$sigC[a$m$sigC]);

  if(newArg){
    a$lNew = ret;
  } else {
    a$lOld = ret;
  }

  a
}

lC = function(a, n, newArg){ # host
  a = lC_kernel1(a, n);
  a = lC_kernel2(a, n);
  a = lC_kernel3(a, n, newArg);
  a
}

lEps = function(a, n, g, arg){

  ret = a$y[n, g] * arg - 
        exp(a$c[a$m$c, n] + arg + mu(a, n, a$phi[a$m$phi, g], 
                                     a$alp[a$m$alp, g], a$del[a$m$del, g])) 
          - arg^2 / (2 * a$eta[a$m$eta, g]^2);

  return(ret);
}

lD = function(a, arg){
  
  if(arg <= 0 || arg > a$d0)
    return(-Inf);
  
  for(g in 1:a$G){ # PARALLELIZE
    a$tmp1[g] = 2 * log(a$eta[a$m$eta, g]);
    a$tmp2[g] = 1/(a$eta[a$m$eta, g] * a$eta[a$m$eta, g]);
  }

  s1 = 0;
  for(g in 1:a$G) # PARALLELIZE
    s1 = s1 + a$tmp1[g];

  s2 = 0;
  for(g in 1:a$G) # PARALLELIZE
    s2 = s2 + a$tmp2[g];

  a$tmp1[1] = arg * a$tau[a$m$tau]^2 / 2;

  ret = -a$G * lgamma(arg/2) + (a$G * arg / 2) * log(a$tmp1[1]);
  ret = ret  - (arg/2 + 1) * s1 - a$tmp1[1] * s2;

  return(ret);
}

lPhi = function(a, g, arg){
 
  s = 0; 
  for(n in 1:a$N){
    a$tmp[0] = mu(a, n, arg, a$alp[a$m$alp, g], a$del[a$m$del, g]);
    s = s + a$y[n, g] * a$tmp[0] - exp(a$c[a$m$c, n] + 
        a$eps[a$m$eps, n, g] + a$tmp[0]);
  }
 
  ret = s - (arg - a$thePhi[a$m$thePhi])^2 / (2 * a$sigPhi[a$m$sigPhi]^2);
  return(ret);
}

lAlp = function(a, g, arg){
  
  s = 0; 
  for(n in 1:a$N){
    if(grp[n] != 2){
      a$tmp[1] = mu(a, n, a$phi[a$m$phi, g], arg, a$del[a$m$del, g]);
      s = s + y[n, g] * a$tmp1[1] - exp(a$c[a$m$c, n] + 
          a$eps[a$m$eps, n, g] + a$tmp1[1]);
    }
  }
 
  if(arg != 0){
    a$tmp1[1] = -(arg - a$theAlp[a$m$theAlp])^2 / (2 * a$sigAlp[a$m$sigAlp]^2) -
                log(1 - a$piAlp[a$m$piAlp]);
  } else {
    a$tmp1[1] = log(a$piAlp[a$m$piAlp])
  }

  ret = s + a$tmp1[1];
  return(ret);
}

lDel = function(a, g, arg){
  
  s = 0; 
  for(n in 1:a$N){
    if(grp[n] != 2){
      a$tmp[1] = mu(a, n, a$phi[a$m$phi, g], a$alp[a$m$alp, g], arg);
      s = s + y[n, g] * a$tmp1[1] - exp(a$c[a$m$c, n] + 
          a$eps[a$m$eps, n, g] + a$tmp1[1]);
    }
  }
 
  if(arg != 0){
    a$tmp1[1] = -(arg - a$theDel[a$m$theDel])^2 / (2 * a$sigDel[a$m$sigDel]^2) -
                log(1 - a$piDel[a$m$piDel]);
  } else {
    a$tmp1[1] = log(a$piDel[a$m$piDel])
  }

  ret = s + a$tmp1[1];
  return(ret);
}



# samplers

sampleC_kernel1 = function(a, n){ # kernel <<<1, 1>>>
  a$old = a$c[a$m$c, n];
  a$new = sampleNormal(a$old, a$tun$c[n]);

  a
}

sampleC_kernel2 = function(a, n){ # kernel <<<1, 1>>>
  lp = min(0, a$lNew - a$lOld);
  lu = log(runif(1));
    
  if(lu < lp){ # accept
    a$c[a$m$c + 1, n] = a$new;
    a$tun$c[n] = a$tun$c[n] * 1.1; # Increase the proposal variance to avoid getting 
                                   # stuck in a mode
  } else { # reject
    a$c[a$m$c + 1, n] = a$old;
    a$tun$c[n] = a$tun$c[n] / 1.1; # If you're rejecting too often, decrease the proposal 
                                   # variance to sample closer to the last accepted value.
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

sampleSigC = function(a){

  for(n in 1:a$N) # PARALLELIZE: 1 BLOCK, N THREADS (LOCKSTEP IS FINE)
    a$tmp1[n] = a$c[a$m$c, n]^2;

  rate = 0;
  for(n in 1:a$N) # PARALLELIZE: PARALLEL SUM IN THRUST
    rate = rate + a$tmp1[n];
  
  shape = (a$N - 1) / 2; 
  rate = rate / 2;
  lb = 1 / a$sigC0^2  

  a$sigC[a$m$sigC + 1] = 1/sqrt(sampleGamma(shape, rate, lb));
  a$m$sigC = a$m$sigC + 1;
  return(a)
}

sampleEps = function(a){ # PARALLELIZE: N BLOCKS, G THREADS PER BLOCK (OR SOMETHING BETTER)
  for(g in 1:a$G){
    for(n in 1:a$N){ 
      old = a$eps[a$m$eps, n, g];
      new = sampleNormal(old, a$tun$eps[n, g]);

      lp = min(0, lEps(a, n, g, new) - lEps(a, n, g, old));
      lu = log(runif(1));
      
      if(lu < lp){ # accept
        a$eps[a$m$eps + 1, n, g] = new;
        a$tun$eps[n, g] = a$tun$eps[n, g] * 1.1; 
      } else { # reject
        a$eps[a$m$eps + 1, n, g] = old;
        a$tun$eps[n, g] = a$tun$eps[n, g] / 1.1;
      }
    }
  }

  a$m$eps = a$m$eps + 1;
  return(a)
}

sampleEta = function(a){

  shape = (a$N + a$d[a$m$d]) / 2; 

  for(g in 1:a$G){ # PARALLELIZE: FIGURE OUT HOW MANY BLOCKS AND THREADS PER BLOCK

    for(n in 1:a$N) # MAYBE PARALLELIZABLE, MAYBE NOT
      a$tmp1[n] = a$eps[a$m$eps, n, g]^2;

    rate = 0;
    for(n in 1:a$N) # MAYBE PARALLELIZE, MAYBE NOT
      rate = rate + a$tmp1[n];
  
    rate = (rate + a$d[a$m$d] * a$tau[a$m$tau] * a$tau[a$m$tau]) / 2; 

    a$eta[a$m$eta + 1, g] = 1/sqrt(sampleGamma(shape, rate));
  }

  a$m$eta = a$m$eta + 1;
  return(a)
}

sampleD = function(a){ 

  old = a$d[a$m$d];
  new = sampleNormal(old, a$tun$d);

  lp = min(0, lD(a, new) - lD(a, old));
  lu = log(runif(1));
    
  if(lu < lp){ # accept
    a$d[a$m$d + 1] = new;
    a$tun$d = a$tun$d * 1.1; # Increase the proposal variance to avoid getting 
                                 # stuck in a mode
  } else { # reject
    a$d[a$m$d + 1] = old;
    a$tun$d = a$tun$d / 1.1; # If you're rejecting too often, decrease the proposal 
                                 # variance to sample closer to the last accepted value.
  }

  a$m$d = a$m$d + 1;
  return(a)
}

sampleTau = function(a){
  shape = a$aTau + a$G * a$d[a$m$d] / 2;

  for(g in 1:a$G) # PARALLELIZE
    a$tmp1[g] = 1/a$eta[a$m$eta, g]^2;

  rate = 0;
  for(g in 1:a$G) # PARALLELIZE
    rate = rate + a$tmp1[g];
  rate = rate * a$d[a$m$d] / 2 + a$bTau;

  a$tau[a$m$tau + 1] = 1/sqrt(sampleGamma(shape, rate));
  a$m$tau = a$m$tau + 1;
  return(a)
}

samplePhi = function(a){ 
  for(g in 1:a$G){ # PARALLELIZE

    old = a$phi[a$m$phi, g];
    new = sampleNormal(old, a$tun$phi[g]);

    lp = min(0, lPhi(a, g, new) - lPhi(a, g, old));
    lu = log(runif(1));
    
    if(lu < lp){ # accept
      a$phi[a$m$phi + 1, g] = new;
      a$tun$phi[g] = a$tun$phi[g] * 1.1; 

    } else { # reject
      a$phi[a$m$phi + 1, g] = old;
      a$tun$phi[g] = a$tun$phi[g] / 1.1; 
    }
  }

  a$m$phi = a$m$phi + 1;
  return(a)
}

sampleThePhi = function(a){
  sm = 0; 
  for(g in 1:a$G) # PARALLELIZE
    sm = sm + a$phi[a$m$phi, g];

  gs = a$gamPhi^2;
  ss = a$sigPhi[a$m$sigPhi]^2;
  den = (a$G * gs + ss);

  m = gs * sm / den;
  s = gs * ss / den;

  a$thePhi[a$m$thePhi + 1] = sampleNormal(m, s);
  a$m$thePhi = a$m$thePhi + 1;
  return(a)
}

sampleSigPhi = function(a){
  shape = (a$G - 1) / 2;

  for(g in 1:a$G) # PARALLELIZE
    a$tmp1[g] = (a$phi[a$m$phi, g] - a$thePhi[a$m$thePhi])^2;

  rate = 0;
  for(g in 1:a$G) # PARALLELIZE
    rate = rate + a$tmp1[g];
  rate = rate / 2;

  lb = 1/a$sigPhi0^2;
  a$sigPhi[a$m$sigPhi + 1] = 1/sqrt(sampleGamma(shape, rate, lb));
  a$m$sigPhi = a$m$sigPhi + 1;
  return(a)
}

sampleAlp = function(a){ 
  for(g in 1:a$G){ # PARALLELIZE

    old = a$alp[a$m$alp, g];

    u = runif(1);
    if(u < a$piAlp[a$m$piAlp]) {
      new = 0;
    } else {

      tmp = 0;
      Nalp = 0;
      for(n in 1:a$N)
        if(a$grp[n] != 2){
          tmp = tmp + a$y[n, g];
          Nalp = Nalp + 1;
        }
      
      mu = a$theAlp[a$m$theAlp] / a$gamAlp^2 + 
           tmp / (Nalp * a$sigAlp[a$m$sigAlp]^2)
      mu = mu / (1/a$gamAlp^2 + 1/a$sigAlp[a$m$sigAlp]^2)

      s = 1/sqrt(1/a$gamAlp^2 + 1/a$sigAlp[a$m$sigAlp]^2);
      new = sampleNormal(mu, s);
    }

    lp = min(0, lAlp(a, g, new) - lAlp(a, g, old));
    lu = log(runif(1));
    
    if(lu < lp){ # accept
      a$alp[a$m$alp + 1, g] = new;
    } else { # reject
      a$alp[a$m$alp + 1, g] = old;
    }
  }

  a$m$alp = a$m$alp + 1;
  return(a)
}

sampleTheAlp = function(a){
  sm = 0;
  Galp = 0;

  for(g in 1:a$G){ # PARALLELIZE
    if(a$alp[a$m$alp, g]){
      a$tmp1[g] = 1;
      a$tmp2[g] = a$alp[a$m$alp, g];
    } else {
      a$tmp1[g] = 0;
      a$tmp2[g] = 0;
    }
  }

  for(g in 1:a$G) # PARALLELIZE
    Galp = Galp + a$tmp1[g];

  for(g in 1:a$G) # PARALLELIZE
    sm = sm + a$tmp2[g];

  gs = a$gamAlp^2;
  ss = a$sigAlp[a$m$sigAlp]^2;
  den = (Galp * gs + ss);

  m = gs * sm / den;
  s = gs * ss / den;

  a$theAlp[a$m$theAlp + 1] = sampleNormal(m, s);
  a$m$theAlp = a$m$theAlp + 1;
  return(a);
}

sampleSigAlp = function(a){
  Galp = 0;
  rate = 0;

  for(g in 1:a$G){ # PARALLELIZE
    if(a$alp[a$m$alp, g]){
      a$tmp1[g] = (a$alp[a$m$alp, g] - a$theAlp[a$m$theAlp])^2;
      a$tmp2[g] = 1;
    } else {
      a$tmp1[g] = 0;
      a$tmp2[g] = 0;
    }
  }

  for(g in 1:a$G) # PARALLELIZE
    rate = rate + a$tmp1[g];

  for(g in 1:a$G) # PARALLELIZE
    Galp = Galp + a$tmp2[g];   

  shape = (Galp - 1) / 2;
  rate = rate / 2;
  lb = 1/a$sigAlp0^2;

  a$sigAlp[a$m$sigAlp + 1] = 1/sqrt(sampleGamma(shape, rate, lb));
  a$m$sigAlp = a$m$sigAlp + 1;
  return(a)
}

samplePiAlp = function(a){ # bookmark
  for(g in 1:a$G){ # PARALLELIZE
    if(a$alp[a$m$alp, g]){
      a$tmp1[g] = 1;
    } else {
      a$tmp1[g] = 0;
    }
  }

  Galp = 0;
  for(g in 1:a$G) # PARALLELIZE
    Galp = Galp + a$tmp1[g];

  a$piAlp[a$m$piAlp + 1] = sampleBeta(a$G + Galp + a$aTau, Galp + a$bTau);
  a$m$piAlp = a$m$piAlp + 1;
  return(a)
}

sampleDel = function(a){ 
  for(g in 1:a$G){ # PARALLELIZE

    old = a$del[a$m$del, g];

    u = runif(1);
    if(u < a$piDel[a$m$piDel]) {
      new = 0;
    } else {

      tmp = 0;
      Ndel = 0;
      for(n in 1:a$N)
        if(a$grp[n] != 2){
          tmp = tmp + a$y[n, g];
          Ndel = Ndel + 1;
        }
      
      mu = a$theDel[a$m$theDel] / a$gamDel^2 + 
           tmp / (Ndel * a$sigDel[a$m$sigDel]^2)
      mu = mu / (1/a$gamDel^2 + 1/a$sigDel[a$m$sigDel]^2)

      s = 1/sqrt(1/a$gamDel^2 + 1/a$sigDel[a$m$sigDel]^2);
      new = sampleNormal(mu, s);
    }

    lp = min(0, lDel(a, g, new) - lDel(a, g, old));
    lu = log(runif(1));
    
    if(lu < lp){ # accept
      a$del[a$m$del + 1, g] = new;
    } else { # reject
      a$del[a$m$del + 1, g] = old;
    }
  }

  a$m$del = a$m$del + 1;
  return(a)
}

sampleTheDel = function(a){

  for(g in 1:a$G){ # PARALLELIZE
    if(a$del[a$m$del, g]){
      a$tmp1[g] = 1;
      a$tmp2[g] = a$del[a$m$del, g];
    } else {
      a$tmp1[g] = 0;
      a$tmp2[g] = 0;
    }
  }

  Gdel = 0;
  for(g in 1:a$G) # PARALLELIZE  
    Gdel = Gdel + a$tmp1[g];

  sm = 0;
  for(g in 1:a$G) # PARALLELIZE
    sm = sm + a$tmp2[g];

  gs = a$gamDel^2;
  ss = a$sigDel[a$m$sigDel]^2;
  den = (Gdel * gs + ss);

  m = gs * sm / den;
  s = gs * ss / den;

  a$theDel[a$m$theDel + 1] = sampleNormal(m, s);
  a$m$theDel = a$m$theDel + 1;
  return(a)
}

sampleSigDel = function(a){
  Gdel = 0;
  rate = 0;

  for(g in 1:a$G){ # PARALLELIZE
    if(a$del[a$m$del, g]){
      a$tmp1[g] = (a$del[a$m$del, g] - a$theDel[a$m$theDel])^2;
      a$tmp2[g] = 1;
    } else {
      a$tmp1[g] = 0;
      a$tmp2[g] = 0;
    }
  }

  for(g in 1:a$G) # PARALLELIZE
    rate = rate + a$tmp1[g];
 
  for(g in 1:a$G) # PARALLELIZE
    Gdel = Gdel + a$tmp2[g];

  shape = (Gdel - 1) / 2;
  rate = rate / 2;
  lb = 1/a$sigDel0^2;

  a$sigDel[a$m$sigDel + 1] = 1/sqrt(sampleGamma(shape, rate, lb));
  a$m$sigDel = a$m$sigDel + 1;
  return(a)
}

samplePiDel = function(a){

  for(g in 1:a$G){ # PARALLELIZE
    if(a$del[a$m$del, g]){
      a$tmp1[g] = 1; 
    } else {
      a$tmp1[g] = 0;
    }
  } 

  Gdel = 0;
  for(g in 1:a$G) # PARALLELIZE
     Gdel = Gdel + a$tmp1[g];

  a$piDel[a$m$piDel + 1] = sampleBeta(a$G + Gdel + a$aTau, Gdel + a$bTau);
  a$m$piDel = a$m$piDel + 1;
  return(a)
}

runChain = function(a){
  for(m in 1:(a$M - 1)){

    print(paste(m))
    print("  step 1")

    a = sampleC(a);

print("  step 2")

    a = sampleTau(a);
    a = samplePiAlp(a);
    a = samplePiDel(a);

print("  step 3")

    a = sampleD(a);
    a = sampleThePhi(a);
    a = sampleTheAlp(a);
    a = sampleTheDel(a);

print("  step 4")

    a = sampleSigC(a); print("    sigc")
    a = sampleSigPhi(a);  print("    sigphi")
    a = sampleSigAlp(a);  print("    sigalp")
    a = sampleSigDel(a);  print("    sigdel")
    a = sampleEta(a); print("    eta")

print("  step 5")

    a = sampleEps(a);

print("  step 6")

    a = samplePhi(a);

print("  step 7")

    a = sampleAlp(a);

print("  step 8")

    a = sampleDel(a);
  }

  return(a)
}

run = function(){
  h = hammer();
  y = h$y
  grp = h$grp
  a = newChain(y, grp, 5, 4, 10)
  return(runChain(a)) 
}