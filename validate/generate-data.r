mu = function(n, phi, alp, del){
  if(n == 1){
    return(phi - alp)
  }
  else if(n == 2){
    return(phi + del)
  } else if(n == 3){
    return(phi + alp)
  } else {
    return(NaN)
  }
}

parms = list(

# initialization constants

sigC0 = 2,

d0 = 100,
aTau = 2,
bTau = 2,

aAlp = 1/2,
aDel = 1/2,

bAlp = 1/2,
bDel = 1/2,

gamPhi = 2,
gamAlp = 2,
gamDel = 2,

sigPhi0 = 2,
sigAlp0 = 2,
sigDel0 = 2
)

parms$G = 100 # number of genes
parms$group = c(1,1,2,2,3,3)
parms$N = length(parms$group)


# iteratively get conditional means/simulated values up the hierarchy
parms$sigC = runif(1, 0, parms$sigC0)
parms$c = rnorm(parms$N, 0, parms$sigC)

parms$d = runif(1, 2, parms$d0)
parms$tau = sqrt(rgamma(1, shape = parms$aTau, rate = parms$bTau))
parms$eta = 1/sqrt(rgamma(parms$G, shape = parms$d/2, rate = parms$d * parms$tau^2/2))
parms$eps = matrix(nrow = parms$G, ncol =parms$ N)
for(g in 1:parms$G)
  for(n in 1:parms$N)
    parms$eps[g, n] = rnorm(1, 0, parms$eta[g])

parms$thePhi = rnorm(1, 0, parms$gamPhi)
parms$sigPhi = runif(1, 0, parms$sigPhi0)
parms$phi = rnorm(parms$G, parms$thePhi, parms$sigPhi)

parms$theAlp = rnorm(1, 0, parms$gamAlp)
parms$sigAlp = runif(1, 0, parms$sigAlp0)
parms$piAlp = rbeta(1, parms$aAlp, parms$bAlp)
parms$alp = rnorm(parms$G, parms$theAlp, parms$sigAlp)

parms$theDel = rnorm(1, 0, parms$gamDel)
parms$sigDel = runif(1, 0, parms$sigDel0)
parms$piDel = rbeta(1, parms$aDel, parms$bDel)
parms$del = rnorm(parms$G, parms$theDel, parms$sigDel)

y = matrix(nrow = parms$G, ncol = parms$N)
for(g in 1:parms$G)
  for(n in 1:parms$N)
    y[g, n] = rpois(1, exp(parms$c[n] + parms$eps[g, n] + mu(parms$group[n], parms$phi[g], parms$alp[g], parms$del[g])))
y = as.data.frame(y)

saveRDS(parms, "parms.rds")
write.table(y, "data.txt", row.names = F, col.names = F)
group = write(parms$group, "group.txt", ncolumns = length(parms$group))
