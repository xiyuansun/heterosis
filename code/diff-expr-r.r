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

emptyChain = function(y, grp, par, M, N, G){
 list(
    yMeanG = rowMeans(y),
    grp = grp,
    M = M,
    N = N,
    G = G,
    
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
      sig = array(0, c(M, G)),
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
    chn$sig[1, g] = 1/sqrt(rgamma(1, shape = chn$d[1] / 2, 
                                  rate = chn$d[1] * chn$tau[1]^2 / 2))


  for(n in 1:N)
    for(g in 1:G)
      chn$eps[1, n, g] = rnorm(1, 0, chn$sig[1, g]);


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

  return(chn);
}

# struct to store the current place in the chain for each parameter

newM = function(){
  list(
    c = 1,
    sigC = 1,

    eps = 1,
    sig = 1,
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

lC = function(y, chn, m, n, arg){
  G = chn$G;

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

lEps = function(y, chn, m, n, g, arg){
  G = chn$G;
  
  y = y[n, g]
  c = chn$c[m$c, n];
  sig = chn$sig[m$sig, g];
  phi = chn$phi[m$phi, g];
  alp = chn$alp[m$alp, g];
  del = chn$del[m$del, g];

  ret = y * arg - exp(c + arg + mu(chn, n, phi, alp, del)) - arg^2 / s * sig^2;
  return(ret);
}

lD = function(y, chn, m, arg){
  G = chn$G;
  d0 = chn$d0;
  tau = chn$tau[m$tau];
  sig = chn$sig[m$sig,];

  if(arg < 0 || arg > d0)
    return(-Inf);

  s1 = 0;
  for(g in 1:G)
    s1 = s1 + 2 * log(sig[g]);

  s2 = 0;
  for(g in 1:G)
    s2 = s2 + (1/sig[g])^2

  tmp = arg * tau^2 / 2

  ret = -G * lgamma(arg/2) + (G * arg / 2) * log(tmp);
  ret = ret  - (arg/2 + 1) * s1 - tmp * s2;

  return(ret);
}

lPhi = function(y, chn, m, g, arg){
  N = chn$N;
  y = y[n, g];

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

lAlp = function(y, chn, m, g, arg){
  N = chn$N;
  y = y[n, g];
  grp = chn$grp

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

lDel = function(y, chn, m, g, arg){
  N = chn$N;
  y = y[n, g];
  grp = chn$grp

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
    tmp = -(arg - theDel)^2 / 2 * sigDel^2 - log(1 - piDel)
  } else {
    tmp = log(piDel)
  }

  ret = s + tmp;
  return(ret);
}


# sampling from known distributions

sampleNormal = function(m = 0, s = 1){
  u1 = runif(1);
  u2 = runif(1);
  return(sqrt(-2 * log(u1)) * sin(2 * 3.14159265 * u2) * s + m);
}

sampleGamma = function(shape = 1, rate = 1, lb = -1){
  if(shape >= 1){ # Marsaglia and Tsang (2000)

    d = shape - 1/3;
    c = 1 / sqrt(9 * d);

    while(1){
      v = -1;
      while(v <= 0){
        x = sampleNormal();
        v = (1 + c*x)^3;
      }

      ret = d * v / rate
      u = runif(1);

      if(u < 1 - 0.0331 * x^4)
        return(ret);

      if(log(u) < 0.5 * x^2 + d * (1 - v + log(v)))
        return(ret);

    }
  } else if (0.1 <= shape && shape < 1){ # Kundu and Gupta (2006)
    
    while(1){
      u = runif(1);
      x = -2 * log(1 - u^(1 / shape));
      v = runif(1);

      tmp1 = exp(-x/2);
      tmp2 = x^(shape - 1)* tmp1 * 2^(1 - shape) * (1 - tmp1)^(1 - shape);

      if(x != 0)
        if(v < tmp2)
          return(x / rate);

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
        z = log(runif(1))/lam
      }
      
      if(z >= 0){
        haznaz = exp(-exp(-z / shape));
      } else{
        haznaz = 1/(w * lam) * exp((lam - 1) * z -exp(z / shape));
      }

      if(haznaz > runif(1))
        return(exp(-z/shape) / rate);
    }
  }
}

shape = .145
pvs = c()
for(j in 1:1000){
  r = c()
  for(i in 1:1000)
    r = c(r, sampleGamma(shape = shape))
  pvs = c(pvs, ks.test(r, f)$p.value)
}
hist(pvs)


trapInt = function(f, a, b, n = 1000){

  sum = 0;
  h = (b - a) / n;

  for(i in 1:n){
    sum = sum + h * f(a + (i + 0.5)*h);
  }

  return(sum);
}

mcInt(h, a, b, n = 1000, fmean, fsd){
  


}


m = 2.5
a = 1
b = 3

e = 0.5
r = c()
for(i in 1:1000){
  u = runif(1)
  if(u < (m - a) / (b - a))
    r = c(r, u^(1 - e) * (b - a) + a)
  else
    r = c(r, u^(1 + e) * (b - a) + a)
}