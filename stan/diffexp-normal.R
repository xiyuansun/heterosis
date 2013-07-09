library(rstan)
library(plyr)


# Simulate data
G = 100
k = rep(c(-1,1), each=3)

# Gene specific parameters
phi = rnorm(G)   # Gene averages
alpha = rnorm(G) # Gene differences
#alpha[1:floor(G/10)] = rnorm(floor(G/10)) # to set some to zero

# Data
d = expand.grid(treatment = k, gene=1:G)
d$y = rpois(nrow(d), exp(phi[d$gene]+alpha[d$gene]*d$treatment))

# Determine genes for which one treatment is always zero
zeros = function(d) {
  prod(by(d$y, d$treatment, sum))==0
}
z = daply(d, "gene", zeros)

## MCMC
diffexp_model = '
  data {
    int<lower=1> G;       // number of genes
    int<lower=1> N;       // number of total observations
    int<lower=0> y[N];    // gene expression
    int<lower=1> gene[N]; // gene ID (1:G)
    real treatment[N];    // treatment ID (-1 and 1)
  }
  parameters {
    real phi[G];
    real alpha[G];
    real mug;
    real mua;
    real<lower=0> sigmag;
    real<lower=0> sigmaa;
  } 
  model {
    for (i in 1:N) {
      y[i] ~ poisson(exp(phi[gene[i]]+alpha[gene[i]]*treatment[i]));
    }

    for (g in 1:G) {
      phi[g]   ~ normal(mug,sigmag);
      alpha[g] ~ normal(mua,sigmaa);
    }

    mug ~ normal(0,1);
    mua ~ normal(0,1);
    sigmag ~ uniform(0,1);
    sigmaa ~ uniform(0,1);
  }
'

diffexp_data = list(G=G, N=nrow(d), y=d$y, gene=d$gene, treatment=d$treatment)

fit <- stan(model_code = diffexp_model, data=diffexp_data, iter=10000, chains=4)
fit1 <- stan(fit = fit, data=diffexp_data, iter=10000, chains=4)

s = monitor(fit1)


# Single gene analysis
onegene = function(d) {
  m = glm(y~as.factor(treatment), poisson, d)
  return(coef(m)[2]/2)
}
est = daply(d, "gene", onegene)


alphas = grep("alpha",rownames(s))

# Determine which are nonzero
sig = s[alphas,4] * s[alphas,8] < 0 

par(mfrow=c(1,2))

plot(s[alphas, 6], alpha, cex=0.5, col=ifelse(sig, "black","red"),
     pch=ifelse(z, 4, 19),
     xlim=range(s[alphas,c(4,8)]),
     ylab="True alpha", xlab="Hierarchically estimated alpha")
abline(0,1, col="gray")
abline(v=0, col="gray")
segments(s[alphas, 4], alpha, s[alphas, 8], alpha, col=ifelse(sig, "black","red"))


plot(s[alphas, 6], est, cex=0.5, col=ifelse(sig, "black","red"),
     pch=ifelse(z, 4, 19),
     xlim=range(est), ylim=range(est),
     ylab="Single-gene analysis", xlab="Hierarchically estimated alpha")
abline(0,1, col="gray")
abline(v=0, col="gray")
segments(s[alphas, 4], est, s[alphas, 8], est, col=ifelse(sig, "black","red"))





