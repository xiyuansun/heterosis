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

hammer = function(){
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

safeLog = function(x){ # for log counts
  if(x > 0) {
    return(log(x));
  } else {
    return(log(0.1));
  }
}

# sampling from known distributions

sampleNormal = function(m = 0, s = 1){
  u1 = runif(1);
  u2 = runif(1);
  return(sqrt(-2 * log(u1)) * sin(2 * 3.14159265 * u2) * s + m);
}

sampleGamma = function(shape = 1, rate = 1, lb = 0){
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

sampleBeta = function(a, b){
  x = sampleGamma(a, 1, 0);
  y = sampleGamma(b, 1, 0);
  return(x / (x + y));
}

# data structures and initialization

initParams = function(){
  list(
    sigC0 = 10,
    d0 = 1e3,

    aTau = 1e2,
    aAlp = 1,
    aDel = 1,

    bTau = 1e2,
    bAlp = 1,
    bDel = 1,
  
    gamPhi = 2,
    gamAlp = 2,
    gamDel = 2,

    sigPhi0 = 2,
    sigAlp0 = 2,
    sigDel0 = 2
  );
}

newMs = function(){
  list(
    c = 1,
    sigC = 1,

    eps = 1,
    eta = 1,
    d = 1,
    tau = 1,

    phi = 1,
    alp = 1,
    del = 1,

    thePhi = 1,
    theAlp = 1,
    theDel = 1,

    sigPhi = 1,
    sigAlp = 1,
    sigDel = 1,

    piAlp = 1,
    piDel = 1
  );
}

emptyChain = function(y, grp, par, M, N, G){
 list(
    y = y,
    yMeanG = rowMeans(y),
    grp = grp,
    M = M,
    N = N,
    G = G,

    ms = newMs(),
    
    sigC0 = par$sigC0,
    d0 = par$d0,

    aTau = par$aTau,
    aAlp = par$aAlp,
    aDel = par$aDel,

    bTau = par$bTau,
    bAlp = par$bAlp,
    bDel = par$bDel,
  
    gamPhi = par$gamPhi,
    gamAlp = par$gamAlp,
    gamDel = par$gamDel,

    sigPhi0 = par$sigPhi0,
    sigAlp0 = par$sigAlp0,
    sigDel0 = par$sigDel0,

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
      piDel = rep(0, M)
  );
}

newChain = function(y, grp, M, N, G){
  par = initParams();
  chn = emptyChain(y, grp, par, M, N, G);

  # compute initial values, mostly using priors

  for(n in 1:N)
    chn$c[1, n] = log(quantile(y[n,], .75))
  chn$c[1, ] = chn$c[1, ] - mean(chn$c[1, ])

  chn$sigC[1] = runif(1, 0, chn$sigC0)

  chn$d[1] = runif(1, 0, chn$d0)
  chn$tau[1] = sqrt(rgamma(1, shape = chn$aTau, rate = chn$bTau))

  for(g in 1:G)
    chn$eta[1, g] = 1/sqrt(rgamma(1, shape = chn$d[1] / 2, 
                                  rate = chn$d[1] * chn$tau[1]^2 / 2))


  for(n in 1:N)
    for(g in 1:G)
      chn$eps[1, n, g] = rnorm(1, 0, chn$eta[1, g]);


  chn$thePhi[1] = rnorm(1, 0, chn$gamPhi)
  chn$theAlp[1] = rnorm(1, 0, chn$gamAlp)
  chn$theDel[1] = rnorm(1, 0, chn$gamDel)

  chn$sigPhi[1] = runif(1, 0, chn$sigPhi0)
  chn$sigAlp[1] = runif(1, 0, chn$sigAlp0)
  chn$sigDel[1] = runif(1, 0, chn$sigDel0)

  chn$piAlp[1] = rbeta(1, chn$aAlp, chn$bAlp)
  chn$piDel[1] = rbeta(1, chn$aDel, chn$bDel)

  for(g in 1:G){
    chn$phi[1, g] = rnorm(1, chn$thePhi[1], chn$sigPhi[1]);

    u = runif(1);
    if(u < chn$piAlp[1]){
      chn$alp[1, g] = 0;
    } else {
      chn$alp[1, g] = rnorm(1, chn$theAlp[1], chn$sigAlp[1]);
    }
    
    u = runif(1);
    if(u < chn$piDel[1]){
      chn$del[1, g] = 0;
    } else {
      chn$del[1, g] = rnorm(1, chn$theDel[1], chn$sigDel[1]);
    }
  }

  # tuning parameters for metropolis steps (std deviations of normal distributions)

  chn$tunC = rep(0, N);
  s = 2 * sd(chn$c[1, ]);
  for(n in 1:N)
    chn$tunC[n] = s;

  chn$tunEps = matrix(0, nrow = N, ncol = G);
  for(g in 1:G){
    s = 2 * sd(chn$eps[1,,g]);
   
    for(n in 1:N)
      chn$tunEps[n, g] = s;
  }

  chn$tunD = chn$d[1];
  
  chn$tunPhi = rep(0, G);
  s = 2* sd(chn$phi[1,]);
  for(g in 1:G)
    chn$tunPhi[g] = s;

  return(chn);
}

# struct to store the current place in the chain for each parameter

mu = function(chn, n, phi, alp, del){
  if(chn$grp[n] == 1){
    return(phi - alp);
  } else if(chn$grp[n] == 2){
    return(phi + del);
  } else if(chn$grp[n] == 3){
    return(phi + alp);
  }
}

# log full conditionals with no convenient form

lC = function(chn, n, arg){
  G = chn$G;
  m = chn$ms;

  ybar = chn$yMeanG;
  sigC = chn$sigC[m$sigC];
  eps = chn$eps[m$eps,n,];
  phi = chn$phi[m$phi,];
  alp = chn$alp[m$alp,];
  del = chn$del[m$del,];

  s = 0;
  for(g in 1:G)
    s = s + exp(eps[g] + mu(chn, n, phi[g], alp[g], del[g]));

  ret = arg * G * ybar[n] - exp(arg) * s - arg^2 / 2 * sigC^2;
  return(ret);
}

lEps = function(chn, n, g, arg){
  G = chn$G;
  m = chn$ms;
  
  y = chn$y[n, g]
  c = chn$c[m$c, n];
  eta = chn$eta[m$eta, g];
  phi = chn$phi[m$phi, g];
  alp = chn$alp[m$alp, g];
  del = chn$del[m$del, g];

  ret = y * arg - exp(c + arg + mu(chn, n, phi, alp, del)) - arg^2 / s * eta^2;
  return(ret);
}

lD = function(chn, arg){
  m = chn$ms;
  G = chn$G;
  d0 = chn$d0;
  tau = chn$tau[m$tau];
  eta = chn$eta[m$eta,];

  if(arg < 0 || arg > d0)
    return(-Inf);

  tmp = rep(0, G)
  for(g in 1:G) # PARALLELIZE
    tmp[g] = log(eta[g]);

  s1 = 0;
  for(g in 1:G) # PARALLELIZE
    s1 = s1 + 2 * tmp[g];

  tmp = rep(0, G);
  for(g in 1:G) # PARALLELIZE
    tmp[g] = 1/(eta[g] * eta[g]);

  s2 = 0;
  for(g in 1:G) # PARALLELIZE
    s2 = s2 + tmp[g];

  tmp = arg * tau^2 / 2

  ret = -G * lgamma(arg/2) + (G * arg / 2) * log(tmp);
  ret = ret  - (arg/2 + 1) * s1 - tmp * s2;

  return(ret);
}

lPhi = function(chn, g, arg){
  y = chn$y[n, g];

  N = chn$N;
  m = chn$ms;

  c = chn$c[m$c, n];
  alp = chn$alp[m$alp, g];
  del = chn$del[m$del, g];

  thePhi = chn$thePhi[m$thePhi];
  sigPhi = chn$sigPhi[m$sigPhi];
  
  s = 0; 
  for(n in 1:N){
    tmp = mu(chn, n, arg, alp, del);
    s = s + y * tmp - exp(c + eps + tmp);
  }
 
  ret = s - (arg - thePhi)^2 / 2 * sigPhi^2;
  return(ret);
}

lAlp = function(chn, g, arg){
  y = chn$y[n, g];
  grp = chn$grp;

  m = chn$ms;
  N = chn$N;

  c = chn$c[m$c, n];
  phi = chn$phi[m$phi, g];
  del = chn$del[m$del, g];

  theAlp = chn$theAlp[m$theAlp];
  sigAlp = chn$sigAlp[m$sigAlp];

  piAlp = chn$piAlp[m$piAlp];
  
  s = 0; 
  for(n in 1:N){
    if(grp[n] != 2){
      tmp = mu(chn, n, phi, arg, del);
      s = s + y * tmp - exp(c + eps + tmp);
    }
  }
 
  if(arg != 0){
    tmp = -(arg - theAlp)^2 / 2 * sigAlp^2 - log(1 - piAlp)
  } else {
    tmp = log(piAlp)
  }

  ret = s + tmp;
  return(ret);
}

lDel = function(chn, g, arg){
  y = chn$y[n, g];
  grp = chn$grp;

  N = chn$N;
  m = chn$ms;

  c = chn$c[m$c, n];
  phi = chn$phi[m$phi, g];
  alp = chn$alp[m$alp, g];

  theDel = chn$theDel[m$theDel];
  sigDel = chn$sigDel[m$sigDel];

  piDel = chn$piDel[m$piDel];
  
  s = 0; 
  for(n in 1:N){
    if(grp[n] == 2){
      tmp = mu(chn, n, phi, alp, arg);
      s = s + y * tmp - exp(c + eps + tmp);
    }
  }
 
  if(arg != 0){
    tmp = -(arg - theDel)^2 / 2 * sigDel^2 - log(1 - piDel);
  } else {
    tmp = log(piDel);
  }

  ret = s + tmp;
  return(ret);
}

# samplers

sampleC = function(chn){ # PARALLELIZE: N BLOCKS, 1 THREAD PER BLOCK (NOTHING IN LOCKSTEP)
  for(n in 1:chn$N){ 
    old = chn$c[chn$ms$c, n];
    new = sampleNormal(old, chn$tunC[n]);

    lp = min(0, lC(chn, n, new) - lC(chn, n, old));
    lu = log(runif(1));
    
    if(lu < lp){ # accept
      chn$c[chn$ms$c + 1, n] = new;
      chn$tunC = chn$tunC * 1.1; # Increase the proposal variance to avoid getting 
                                   # stuck in a mode
    } else { # reject
      chn$c[chn$ms$c + 1, n] = new;
      chn$tunC = chn$tunC / 1.1; # If you're rejecting too often, decrease the proposal 
                                   # variance to sample closer to the last accepted value.
    }
  }

  chn$ms$c = chn$ms$c + 1;
}

sampleSigC = function(chn){
  N = chn$N;

  shape = (N - 1) / 2;
  tmp = rep(0, N);  

  for(n in 1:N) # PARALLELIZE: 1 BLOCK, N THREADS (LOCKSTEP IS FINE)
    tmp = chn$c[chn$ms$c, n]^2;

  rate = 0;
  for(n in 1:N){ # PARALLELIZE: PARALLEL SUM IN THRUST
    rate = rate + tmp[n];
  }
  rate = rate / 2;

  lb = 1 / chn$sigC0^2  

  chn$sigC[chn$ms$sigC + 1] = 1/sqrt(sampleGamma(shape, rate, lb));
  chn$ms$sigC = chn$ms$sigC + 1;
}

sampleEps = function(chn){ # PARALLELIZE: N BLOCKS, G THREADS PER BLOCK (OR SOMETHING BETTER)
  for(g in 1:chn$G){
    for(n in 1:chn$N){ 
      old = chn$eps[chn$ms$eps, n, g];
      new = sampleNormal(old, chn$tunEps[n, g]);

      lp = min(0, lEps(chn, n, g, new) - lEps(chn, n, g, old));
      lu = log(runif(1));
      
      if(lu < lp){ # accept
        chn$Eps[chn$ms$Eps + 1, n, g] = new;
        chn$tuneEps[n, g] = chn$tuneEps[n, g] * 1.1; 
      } else { # reject
        chn$Eps[chn$ms$Eps + 1, n] = new;
        chn$tuneEps[n, g] = chn$tuneEps[n, g] / 1.1;
      }
    }
  }

  chn$ms$Eps = chn$ms$Eps + 1;
}

sampleEta = function(chn){
  for(g in 1:chn$G){ # PARALLELIZE: FIGURE OUT HOW MANY BLOCKS AND THREADS PER BLOCK

    shape = (chn$N - chn$d[chn$ms$d]) / 2;
    tmp = rep(0, chn$N);  

    for(n in 1:chn$N) # MAYBE PARALLELIZABLE, MAYBE NOT
      tmp = chn$eps[chn$ms$eps, n, g]^2;

    rate = 0;
    for(n in 1:chn$N) # MAYBE PARALLELIZE, MAYBE NOT
      rate = rate + tmp[n];
  
    rate = (rate + chn$d[chn$ms$d] * chn$tau[chn$ms$tau] * chn$tau[chn$ms$tau]) / 2; 

    chn$eta[chn$ms$eta + 1, g] = 1/sqrt(sampleGamma(shape, rate));
  }

  chn$ms$eta = chn$ms$eta + 1;
}

sampleD = function(chn){ 

  old = chn$d[chn$ms$d];
  new = sampleNormal(old, chn$tunD);

  lp = min(0, lD(chn, new) - lD(chn, old));
  lu = log(runif(1));
    
  if(lu < lp){ # accept
    chn$d[chn$ms$d + 1] = new;
    chn$tunD = chn$tunD * 1.1; # Increase the proposal variance to avoid getting 
                                 # stuck in a mode
  } else { # reject
    chn$d[chn$ms$d + 1] = new;
    chn$tunD = chn$tunD / 1.1; # If you're rejecting too often, decrease the proposal 
                                 # variance to sample closer to the last accepted value.
  }

  chn$ms$d = chn$ms$d + 1;
}

sampleTau = function(chn){
  shape = chn$aTau + chn$G * chn$d[chn$ms$d] / 2;

  tmp = rep(0, G);
  for(g in 1:G) # PARALLELIZE
    tmp[g] = 1/chn$eta[chn$ms$eta, g]^2;

  rate = 0;
  for(g in 1:G) # PARALLELIZE
    rate = rate + tmp[g];
  rate = rate * chn$d[chn$ms$d] / 2 + chn$bTau;

  chn$tau[chn$ms$tau + 1] = 1/sqrt(sampleGamma(shape, rate));
  chn$ms$tau = chn$ms$tau + 1;
}

samplePhi = function(chn){ 
  for(g in 1:chn$G){ # PARALLELIZE

    old = chn$phi[chn$ms$phi, g];
    new = sampleNormal(old, chn$tunPhi[g]);

    lp = min(0, lPhi(chn, g, new) - lPhi(chn, g, old));
    lu = log(runif(1));
    
    if(lu < lp){ # accept
      chn$phi[chn$ms$phi + 1, g] = new;
      chn$tunPhi[g] = chn$tunPhi[g] * 1.1; 

    } else { # reject
      chn$phi[chn$ms$phi + 1, g] = new;
      chn$tunPhi[g] = chn$tunPhi[g] / 1.1; 
    }
  }

  chn$ms$phi = chn$ms$phi + 1;
}

sampleThePhi = function(chn){
  sm = 0; 
  for(g in 1:chn$G) # PARALLELIZE
    sm = sm + chn$phi[chn$ms$phi, g];

  gs = chn$gamPhi^2;
  ss = chn$sigPhi^2;
  den = (chn$G * gs + ss);

  m = gs * sm / den;
  s = sg * ss / den;

  chn$thePhi[chn$ms$thePhi + 1] = sampleNormal(m, s);
  chn$ms$thePhi = chn$ms$thePhi + 1;
}

sampleSigPhi = function(chn){
  shape = (chn$G - 1) / 2;

  rate = 0;
  for(g in 1:G) # PARALLELIZE
    rate = rate + (chn$phi[chn$ms$phi, g] - chn$thePhi[chn$ms$thePhi])^2;
  rate = rate / 2;

  lb = 1/chn$sigPhi0^2;
  chn$sigPhi[chn$ms$sigPhi + 1] = 1/sqrt(sampleGamma(shape, rate, lb));
  chn$ms$sigPhi = chn$ms$sigPhi + 1;
}

sampleAlp = function(chn){ 
  for(g in 1:chn$G){ # PARALLELIZE

    old = chn$alp[chn$ms$alp, g];

    u = unif(1);
    if(u < chn$piAlp[chn$ms$piAlp]) {
      new = 0;
    } else {

      tmp = 0;
      Nalp = 0;
      for(n in 1:chn$N){
        if(chn$grp[n] != 2){
          tmp = tmp + chn$y[n, g];
          Nalp = Nalp + 1;
        }
      
      mu = chn$theAlp[chn$ms$theAlp] / chn$gamAlp^2 + 
           tmp / (Nalp * chn$sigAlp[chn$ms$sigAlp]^2)
      mu = mu / (1/chn$gamAlp^2 + 1/chn$sigAlp[chn$ms$sigAlp]^2)

      s = 1/sqrt(1/chn$gamAlp^2 + 1/chn$sigAlp^2);
      new = sampleNormal(mu, s);
    }

    lp = min(0, lAlp(chn, g, new) - lAlp(chn, g, old));
    lu = log(runif(1));
    
    if(lu < lp){ # accept
      chn$alp[chn$ms$alp + 1, g] = new;
      chn$tunAlp[g] = chn$tunAlp[g] * 1.1; 

    } else { # reject
      chn$alp[chn$ms$alp + 1, g] = new;
      chn$tunAlp[g] = chn$tunAlp[g] / 1.1; 
    }
  }

  chn$ms$alp = chn$ms$alp + 1;
}

sampleTheAlp = function(chn){
  sm = 0;
  Galp = 0;

  for(g in 1:chn$G) # PARALLELIZE
    if(chn$alp[chn$ms$alp, g]){
      Galp = Galp + 1;
      sm = sm + chn$alp[chn$ms$alp, g];
    }

  gs = chn$gamAlp^2;
  ss = chn$sigAlp^2;
  den = (Galp * gs + ss);

  m = gs * sm / den;
  s = sg * ss / den;

  chn$theAlp[chn$ms$theAlp + 1] = sampleNormal(m, s);
  chn$ms$theAlp = chn$ms$theAlp + 1;
}

sampleSigAlp = function(chn){
  Galp = 0;
  rate = 0;

  for(g in 1:G) # PARALLELIZE
    if(chn$alp[chn$ms$alp, g]){
      rate = rate + (chn$alp[chn$ms$alp, g] - chn$theAlp[chn$ms$theAlp])^2;
      Galp = Galp + 1;
    }

  shape = (Galp - 1) / 2;
  rate = rate / 2;
  lb = 1/chn$sigAlp0^2;

  chn$sigAlp[chn$ms$sigAlp + 1] = 1/sqrt(sampleGamma(shape, rate, lb));
  chn$ms$sigAlp = chn$ms$sigAlp + 1;
}

samplePiAlp = function(chn){
  Galp = 0;
  for(g in 1:G) # PARALLELIZE
    if(chn$alp[chn$ms$alp, g])
      Galp = Galp + 1;

  chn$piAlp[chn$ms$piAlp + 1] = sampleBeta(chn$G + Galp + chn$aTau, Galp + chn$bTau);
  chn$ms$piAlp = chn$ms$piAlp + 1;
}

sampleDel = function(chn){ 
  for(g in 1:chn$G){ # PARALLELIZE

    old = chn$del[chn$ms$del, g];

    u = unif(1);
    if(u < chn$piDel[chn$ms$piDel]) {
      new = 0;
    } else {

      tmp = 0;
      Ndel = 0;
      for(n in 1:chn$N){
        if(chn$grp[n] != 2){
          tmp = tmp + chn$y[n, g];
          Ndel = Ndel + 1;
        }
      
      mu = chn$theDel[chn$ms$theDel] / chn$gamDel^2 + 
           tmp / (Ndel * chn$sigDel[chn$ms$sigDel]^2)
      mu = mu / (1/chn$gamDel^2 + 1/chn$sigDel[chn$ms$sigDel]^2)

      s = 1/sqrt(1/chn$gamDel^2 + 1/chn$sigDel^2);
      new = sampleNormal(mu, s);
    }

    lp = min(0, lDel(chn, g, new) - lDel(chn, g, old));
    lu = log(runif(1));
    
    if(lu < lp){ # accept
      chn$del[chn$ms$del + 1, g] = new;
      chn$tunDel[g] = chn$tunDel[g] * 1.1; 

    } else { # reject
      chn$del[chn$ms$del + 1, g] = new;
      chn$tunDel[g] = chn$tunDel[g] / 1.1; 
    }
  }

  chn$ms$del = chn$ms$del + 1;
}

sampleTheDel = function(chn){
  sm = 0;
  Gdel = 0;

  for(g in 1:chn$G) # PARALLELIZE
    if(chn$del[chn$ms$del, g]){
      Gdel = Gdel + 1;
      sm = sm + chn$del[chn$ms$del, g];
    }

  gs = chn$gamDel^2;
  ss = chn$sigDel^2;
  den = (Gdel * gs + ss);

  m = gs * sm / den;
  s = sg * ss / den;

  chn$theDel[chn$ms$theDel + 1] = sampleNormal(m, s);
  chn$ms$theDel = chn$ms$theDel + 1;
}

sampleSigDel = function(chn){
  Gdel = 0;
  rate = 0;

  for(g in 1:G) # PARALLELIZE
    if(chn$del[chn$ms$del, g]){
      rate = rate + (chn$del[chn$ms$del, g] - chn$theDel[chn$ms$theDel])^2;
      Gdel = Gdel + 1;
    }

  shape = (Gdel - 1) / 2;
  rate = rate / 2;
  lb = 1/chn$sigDel0^2;

  chn$sigDel[chn$ms$sigDel + 1] = 1/sqrt(sampleGamma(shape, rate, lb));
  chn$ms$sigDel = chn$ms$sigDel + 1;
}

samplePiDel = function(chn){
  Gdel = 0;
  for(g in 1:G) # PARALLELIZE
    if(chn$del[chn$ms$del, g])
      Gdel = Gdel + 1;

  chn$piDel[chn$ms$piDel + 1] = sampleBeta(chn$G + Gdel + chn$aTau, Gdel + chn$bTau);
  chn$ms$piDel = chn$ms$piDel + 1;
}

runChain = function(chn){
  for(m in 1:(chn$M - 1)){

    sampleC(chn);

    sampleTau(chn);
    samplePiAlp(chn);
    samplePiDel(chn);

    sampleD(chn);
    sampleThePhi(chn);
    sampleTheAlp(chn);
    sampleTheDel(chn);

    sampleSigC(chn);
    sampleSigPhi(chn);
    sampleSigAlp(chn);
    sampleSigDel(chn);
    sampleEta(chn);

    sampleEps(chn);

    samplePhi(chn);

    sampleAlp(chn);

    sampleDel(chn);
  }
}