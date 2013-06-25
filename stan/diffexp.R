
diffexp_timing = function(G) {

# Simulate data
#G = 30000
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

compile_time = system.time(fit <- stan(model_code = diffexp_model, data=diffexp_data, iter=4, chains=4))
fit_time = system.time(fit1 <- stan(fit = fit, data=diffexp_data, iter=4, chains=4))

return(data.frame(compile_time=compile_time[3], fit_time=fit_time[3]))

}

Gs = 3*10^(2:4)
timing = matrix(NA, nrow=length(Gs), ncol=3)

for (i in 1:length(Gs)) {
  tmp = diffexp_timing(Gs[i])
  timing[i,] = c(Gs[i], tmp$compile_time, tmp$fit_time)
}

timing[,2:3] = timing[,2:3]/60

plot(timing[,1], timing[,3], log='xy', ylim=range(timing[,2:3]), lwd=2, type="l",
     xlab='Number of genes', ylab='Time (minutes)', main='Elapsed time for Stan')
lines(timing[,1], timing[,2], col=2, lwd=2)
legend('topleft', inset=0.01, c("Compile","Execution"), col=1:2, lwd=2)





