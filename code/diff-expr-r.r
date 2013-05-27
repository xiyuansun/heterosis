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

safelog = function(x){
  if(x > 0) {
    return(log(x));
  } else {
    return(0);
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
      aux = c(phis, safelog(d$y[n, g]) - c$e[1, n, g] - c$c[1, n]);

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

mu = function(m, n, g, c){
  if(c$gr[n] == 1){
    return(c$ph[m, g] - c$al[m, g]);
  } else if(c$gr[n] == 2){
    return(c$ph[m, g] + c$de[m, g]);
  } else if(c$gr[n] == 3){
    return(c$ph[m, g] + c$al[m, g]);
  }
}

# log full conditionals (up to normalizing constants)

l_c = function(m, n, c, d){ # parallelize across genes
  g = 0;
  l = 0;

  for(g in 1:c$G)
    l = l + safelog(dpois(d$y[n, g], exp(c$c[m, n] + c$e[m, n, g] + mu(m, n, g, c))));

  l = l + c$G * safelog(dnorm(c$c[m, n], 0, c$sc[m]));

  return(l);
}

l_sc = function(m, c, d){
  n = 0;
  l = 0;

  if(c$sc[m] == 0)
    return(0);

  for(n in 1:c$N)
    l = l + safelog(dnorm(c$c[m,n], 0, c$sc[m]))

  return(G * l)
}


# d = hammer()
c = new_chain(d, 5)