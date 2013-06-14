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
library(rjags)

hammer = function(){
  load(url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/hammer_eset.RData"))
#  load("../data/hammer_eset.RData")
  counts = exprs(hammer.eset)
  counts = counts[rowSums(counts) > 0,]
  
  group = as.vector(phenoData(hammer.eset)$protocol)
  group[group == "control"] = "1"
  group[group == "L5 SNL"] = "2"
  group = as.factor(as.numeric(group))
  
  list(counts = counts, group = group)
}

l = hammer()
y = l$counts
group = l$group
G = dim(y)[1]
N = dim(y)[2]
c = apply(y, 2, function(x){quantile(x, 0.75)})
c = c / prod(c)^(1/N) 

model = "
model{
  for(g in 1:G){
    for(n in 1:N){
      y[g,n] ~ dpois(exp(mu[g, group[n]] + epsilon[g,n] + log_c[n]))
      epsilon[g,n] ~ dnorm(0, tau[g])
    }

    tau[g] ~ dgamma(d*tau0/2, d/2)

    mu[g, 1] <- phi[g] - alpha[g]
    mu[g, 2] <- phi[g] + delta[g]
#    mu[g, 3] <- phi[g] + alpha[g] # uncomment for heterosis

    phi[g] ~ dnorm(theta_phi, tau_phi)

    alpha[g] <- (1-z_alpha[g])*alpha0[g]
    z_alpha[g] ~ dbern(pi_alpha)
    alpha0[g] ~ dnorm(theta_alpha, tau_alpha)

    delta[g] <- (1-z_delta[g])*delta0[g]
    z_delta[g] ~ dbern(pi_delta)
    delta0[g] ~ dnorm(theta_delta, tau_delta)
  }

  for(n in 1:N){
    log_c[n] ~ dnorm(0, tau_c)
  }

  d ~ dunif(0,1e3)
  tau0 ~ dgamma(1,1)

  pi_alpha ~ dunif(0,1)
  pi_delta ~ dunif(0,1)

  theta_phi ~ dnorm(0, 1e3)
  theta_alpha ~ dnorm(0, 1e3)
  theta_delta ~ dnorm(0, 1e3)

  tau_phi <- 1/sigma_phi^2
  tau_alpha <- 1/sigma_alpha^2
  tau_delta <- 1/sigma_delta^2
  tau_c <- 1/sigma_c^2

  sigma_phi ~ dunif(0, 1e3)
  sigma_alpha ~ dunif(0, 1e3)
  sigma_delta ~ dunif(0, 1e3)
  sigma_c ~ dunif(0, 1)
}
"

d = list(y = y[1:40,], group = group, G = 40, N = N)
m = jags.model(textConnection(model), d, n.chains = 3, n.adapt = 100)
res = coda.samples(m, "phi", 10000)

gelman.diag(res)
plot(res)

