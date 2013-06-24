
# Simulate data
G = 30000
k = rep(c(-1,1), each=3)

# Gene specific parameters
phi = rnorm(G)   # Gene averages
alpha = rep(0,G) # Gene differences
alpha[1:floor(G/10)] = rnorm(floor(G/10))

# Data
d = expand.grid(gene=1:G, treatment=k)
d$y = rpois(nrow(d), exp(phi[d$gene]+alpha[d$gene]*d$treatment))

## MCMC
diffexp_model = '
  data {
    int<lower=0> G;       // number of genes
    int<lower=0> N;       // number of total observations
    int<lower=0> y[N];    // gene expression
    int<lower=0> gene[N]; // gene ID (1:G)
    real treatment[N];    // treatment ID (-1 and 1)
  }
  parameters {
    real phi[G];
    real alpha[G];
  } 
  model {
    for (i in 1:N) {
      y[i] ~ poisson(exp(phi[gene[i]]+alpha[gene[i]]*treatment[i]));
    }

    for (g in 1:G) {
      phi[g]   ~ normal(0,1);
      alpha[g] ~ normal(0,1);
    }
  }
'

diffexp_data = list(G=G, N=nrow(d), y=d$y, gene=d$gene, treatment=d$treatment)

fig = stan(model_code = diffexp_model, data=diffexp_data, iter=4, chains=4)



