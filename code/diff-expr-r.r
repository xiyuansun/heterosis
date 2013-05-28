# Author: Will Landau
# PIs: Dan Nettleton, Jarad Niemi, Peng Liu
# Iowa State University
# May 2013
# 
# This code attempts to learn about differential expression
# between 2 treatment groups of an RNA-Seq dataset.
# 
# Model (based on section 2.3.2 of the proposal):
#
# y_gij | lambda_gij ~ Poisson(c_ij lambda_gij)
# epsilon_gij | sigma_g^2 ~ N(0, sigma_g^2)
# sigma^2_g | d, sigma0^2 ~ d sigma0^2 inv-chisquare_d
# c_ij ~ N(1, sigma_c^2)
# 
# log(lambda_gij) = mu_gi + epsilon_gij
# mu_g1 = phi_g - alpha_g
# mu_g2 = phi_g + delta_g
# mu_g3 = phi_g + alpha_g
#
# phi_g | theta_phi, sigma_phi^2 ~ N(theta_phi, sigma_phi^2)
# alpha_g | theta_alpha, sigma_alpha^2, pi_alpha ~ 
#   pi_alpha 1(alpha_g == 0) + (1- pi_alpha) N(theta_alpha, sigma_alpha^2)
# delta_g | theta_delta, sigma_delta^2, pi_delta ~ 
#   pi_delta 1(delta_g == 0) + (1- pi_delta) N(theta_delta, sigma_delta^2)
# 
# d ~ unif(0, 1e3)
# pi_alpha ~ unif(0, 1) 
# pi_delta ~ unif(0,1)
# 
# theta_phi ~ N(0, 1e3)
# theta_alpha ~ N(0, 1e3)
# theta_delta ~ N(0, 1e3)
# 
# sigma0 ~ unif(0, 1e3)
# sigma_phi ~ unif(0, 1e3)  
# sigma_alpha ~ unif(0, 1e3)  
# sigma_delta ~ unif(0, 1e3)
# sigma_c ~ unif(0, 1)

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
       gr = group, 
       G = 30, #dim(counts)[1], 
       N = dim(counts)[2]))
}

safelog = function(x){ # for log counts
  if(x > 0) {
    return(log(x));
  } else {
    return(log(0.1));
  }
}

new_chain = function(d, M){
  n = 0;
  g = 0;
  shape = 0;
  rate = 0;
  coin = 0;

  N = d$N;
  G = d$G;

  # create chain object

  c = list(
    gr = d$gr,
    M = M,
    N = N,
    G = G,
    
    c = array(0, c(M, N)),
      sc = rep(0, M),
    
    e = array(0, c(M, N, G)),
      s = array(0, c(M, G)),
        d = rep(0, M),
        s0 = rep(0, M),

    ph = array(0, c(M, G)),
      th_ph = rep(0, M),
      s_ph = rep(0, M),

    al = array(0, c(M, G)),
      th_al = rep(0, M),
      s_al = rep(0, M),
      pi_al = rep(0, M),

    de = array(0, c(M, G)),
      th_de = rep(0, M),
      s_de = rep(0, M),
      pi_de = rep(0, M)
  );

  # compute initial values, mostly using priors

  c$sc[1] = c(runif(1));

  for(n in 1:N)
    c$c[1, n] = log(quantile(d$y[n,], .75))

  c$c[1, ] = c$c[1, ] - mean(c$c[1, ])

  ##

  c$s0[1] = sqrt(rexp(1));
  c$d[1] = runif(1, 0, 1000);

  shape = c$d[1] * c$s0[1]^2 / 2;
  rate = c$d[1] / 2;

  for(g in 1:G)
    c$s[1, g] = rgamma(1, shape = shape, rate = rate);

  for(n in 1:N)
    for(g in 1:G)
      c$e[1, n, g] = rnorm(1, 0, c$s[1,g]);

  ##

  c$s_ph[1] = runif(1, 0, 1000);
  c$th_ph[1] = rnorm(1, 0, sqrt(1000));

  for(g in 1:G){ # change phi to the parental average for heterosis
    aux = c()

    for(n in 1:N)
      aux = c(aux, safelog(d$y[n, g]) - c$e[1, n, g] - c$c[1, n]);

    c$ph[1, g] = median(aux);
  }

  ##

  c$s_al[1] = runif(1, 0, 1000);
  c$th_al[1] = rnorm(1, 0, sqrt(1000));
  c$pi_al[1] = runif(1);

  for(g in 1:G){ # change alpha to half the parental difference for heterosis
    aux = c();

    for(n in 1:N)
      if(n == 1)
        aux = c(aux, c$ph[1, g] + c$e[1, n, g] + c$c[1, n] - safelog(d$y[n, g]));

    c$al[1, g] = median(aux);
  }

  ##

  c$s_de[1] = runif(1, 0, 1000);
  c$th_de[1] = rnorm(1, 0, sqrt(1000));
  c$pi_de[1] = runif(1);

  for(g in 1:G){ # change delta to the parent-child difference for heterosis
    aux = c();

    for(n in 1:N)
      if(n == 2)
        aux = c(aux, - c$ph[1, g] - c$e[1, n, g] - c$c[1, n] + safelog(d$y[n, g]));

    c$de[1, g] = median(aux);
  }

  return(c); 
}

mu = function(n, ph, al, de){
  if(c$gr[n] == 1){
    return(ph - al);
  } else if(c$gr[n] == 2){
    return(ph + de);
  } else if(c$gr[n] == 3){
    return(ph + al);
  }
}

# log full conditionals (up to normalizing constants)

l_c = function(arg, m, n, c, d){ # parallelize across genes
  g = 0;
  l = 0;

  for(g in 1:c$G)
    l = l + safelog(dpois(d$y[n, g], 
                    exp(arg + c$e[m, n, g] + mu(n, c$ph[m, g], c$al[m, g], c$de[m, g]))));

  l = l + c$G * safelog(dnorm(arg, 0, c$sc[m]));

  return(l);
}

for(n in 1:c$N)
  print(paste(l_sc(c$sc[m], m, c, d)))

l_sc = function(arg, m, c, d){
  n = 0;
  l = 0;

  if(arg < 0 || arg > 1)
    return(-Inf);

  for(n in 1:c$N)
    l = l + safelog(dnorm(c$c[m,n], 0, arg));

  return(G * l);
}

l_e = function(arg, m, n, g, c, d){
  l = 0;
  l = l + safelog(dpois(d$y[n, g], 
                 exp(c$c[m, n] + arg + mu(n, c$ph[m, g], c$al[m, g], c$de[m, g]))));
  return(l + safelog(dnorm(arg, 0, c$s[m, g])));
}

l_s = function(arg, m, g, c, d){
  n = 0;
  l = 0;  

  for(n in 1:c$N)
    l = l + safelog(dnorm(c$e[m, n, g], 0, arg));
  
  l = l + c$N * safelog(dgamma(arg, shape = c$d[m] * c$s0[m]^2 / 2, rate = c$d[m] / 2));

}

l_d = function(arg, m, c, d){ # parallelize accross genes
  g = 0;
  l = 0;

  if(arg < 0 || arg > 1000)
    return(-Inf);

  for(g in 1:c$G)
    l = l + safelog(dgamma(c$s[m, g], shape = arg * c$s0[m]^2 / 2, rate = arg / 2));

  return(l * c$N);
}

l_s0_sq = function(arg, m, c, d){
  g = 0;
  l = 0;

  if(arg < 0)
    return(-Inf);

  for(g in 1:c$G)
    l = l + safelog(dgamma(c$s[m, g], shape = c$d[m] * arg / 2, rate = c$d[m] / 2));

  l = l * c$N;
  return(c$N * c$G * safelog(dexp(arg)));
}

l_ph = function(arg, m, g, c, d){
  n = 0;
  l = 0;

  for(n in 1:c$N)
    l = l + safelog(dpois(d$y[n, g], 
                    exp(c$c[m, n] + c$e[m, n, g] + mu(n, arg, c$al[m, g], c$de[m, g]))));
  
  return(l + c$N * safelog(dnorm(arg, c$th_ph[m], c$s_ph[m])));
}

l_th_ph = function(arg, m, c, d){
  g = 0;
  l = 0;

  for(g in 1:c$G)
    l = l + safelog(dnorm(c$ph[m, g], arg, c$s_ph));

  l = l + c$G * safelog(dnorm(arg, 0, sqrt(1000)));

  return(c$N * l);
}

l_s_ph = function(arg, m, c, d){
  g = 0;
  l = 0;

  if(arg < 0 || arg > 1000)
    return(-Inf);

  for(g in 1:c$G)
    l = l + safelog(dnorm(c$ph[m, g], c$th_ph[m], arg));

  return(c$N * l);
}

l_al = function(arg, m, g, c, d){
  n = 0;
  l = 0;

  for(n in 1:c$N)
    if(c$gr[n] != 2)
      l = l + safelog(dpois(d$y[n, g], 
                      exp(c$c[m, n] + c$e[m, n, g] + mu(n, arg, c$al[m, g], c$de[m, g]))));

  if(arg == 0){
    return(l + sum(c$gr != 2) * safelog(c$pi_al[m]));
  } else {
    return(l + sum(c$gr != 2) * (safelog(1 - c$pi_al[m]) + 
                                 safelog(dnorm(arg, c$th_al[m], c$s_al[m]))));
  }
}

l_th_al = function(arg, m, c, d){
  g = 0;
  n = 0;
  l = 0;

  for(n in 1:c$N)
    if(c$gr[n] != 2)
      for(g in 1:c$G)
        l = l + safelog(dnorm(c$al[m, g], arg, c$s_al));

  return(l + c$G * c$N * safelog(dnorm(arg, 0, sqrt(1000))));
}

l_s_al = function(arg, m, c, d){
  g = 0;
  n = 0;
  l = 0;

  if(arg < 0 || arg > 1000)
    return(-Inf);

  for(n in 1:c$N)
    if(c$gr[n] != 2)
      for(g in 1:c$G)
        l = l + safelog(dnorm(c$al[m, g], c$th_al[m], arg));

  return(l);
}

l_pi_al = function(arg, m, c, d){
  n = 0;
  g = 0;
  l = 0;

  if(arg < 0 || arg > 1)
    return(-Inf);

  s = sum(c$gr != 2);
  
  for(g in 1:c$G){
    if(c$al[m, g] == 0) {
      l = l + safelog(arg);
    } else {
      l = l + safelog(1 - arg) + safelog(dnorm(c$al[m, g], c$th_al[m], c$s_al[m]));
    }
  }

  return(s * l);
}

l_de = function(arg, m, g, c, d){
  n = 0;
  l = 0;

  for(n in 1:c$N)
    if(c$gr[n] == 2)
      l = l + safelog(dpois(d$y[n, g], 
                      exp(c$c[m, n] + c$e[m, n, g] + mu(n, arg, c$de[m, g], c$de[m, g]))));

  if(arg == 0){
    return(l + sum(c$gr == 2) * safelog(c$pi_de[m]));
  } else {
    return(l + sum(c$gr == 2) * (safelog(1 - c$pi_de[m]) + 
                                 safelog(dnorm(arg, c$th_de[m], c$s_de[m]))));
  }
}

l_th_de = function(arg, m, c, d){
  g = 0;
  n = 0;
  l = 0;

  for(n in 1:c$N)
    if(c$gr[n] == 2)
      for(g in 1:c$G)
        l = l + safelog(dnorm(c$de[m, g], arg, c$s_de));

  return(l + c$G * c$N * safelog(dnorm(arg, 0, sqrt(1000))));
}

l_s_de = function(arg, m, c, d){
  g = 0;
  n = 0;
  l = 0;

  if(arg < 0 || arg > 1000)
    return(-Inf);

  for(n in 1:c$N)
    if(c$gr[n] == 2)
      for(g in 1:c$G)
        l = l + safelog(dnorm(c$de[m, g], c$th_de[m], arg));

  return(l);
}

l_pi_de = function(arg, m, c, d){
  n = 0;
  g = 0;
  l = 0;

  if(arg < 0 || arg > 1)
    return(-Inf);

  s = sum(c$gr == 2);
  
  for(g in 1:c$G){
    if(c$de[m, g] == 0) {
      l = l + safelog(arg);
    } else {
      l = l + safelog(1 - arg) + safelog(dnorm(c$de[m, g], c$th_de[m], c$s_de[m]));
    }
  }

  return(s * l);
}

## functions to sample from conditionals

## Gibbs sampler function

## run Gibbs sampler

d = hammer()
c = new_chain(d, 5)