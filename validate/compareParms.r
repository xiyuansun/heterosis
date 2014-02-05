trueParms = readRDS("parms.RDS")
hyper1 = read.table("out/hyper-chain1.txt", head = T)
hyper2 = read.table("out/hyper-chain2.txt", head = T)
hyper3 = read.table("out/hyper-chain3.txt", head = T)

for(i in dim(hyper1)[2]:1)
  if(any(hyper1[,i] == "."))
    hyper1 = hyper1[,-i]

for(i in dim(hyper2)[2]:1)
  if(any(hyper2[,i] == "."))
    hyper2 = hyper2[,-i]
  
for(i in dim(hyper3)[2]:1)
  if(any(hyper3[,i] == "."))
    hyper3 = hyper3[,-i]
 
if(is.null(dim(hyper1)))
  hyper1 = data.frame(parm = hyper1)

if(is.null(dim(hyper1)))
  hyper2 = data.frame(parm = hyper2)

if(is.null(dim(hyper1)))
  hyper3 = data.frame(parm = hyper3)

if(!file.exists("fig"))
  dir.create("fig")

#iters = (length(hyper1[1][[1]]) - 1000):length(hyper1[1][[1]])
iters = 1:length(hyper1[1][[1]])

for(n in names(hyper1)){
  pdf(paste("fig/", n, ".pdf", sep=""))
  plot(x = rep(iters, 3), y = c(hyper1[n][[1]][iters], hyper2[n][[1]][iters], hyper3[n][[1]][iters]), lty = 1, ylab = n, xlab = "iteration", col = "white")
  points(x = iters, y = hyper1[n][[1]][iters], col = 1)
  points(x = iters, y = hyper2[n][[1]][iters], col = 2)
  points(x = iters, y = hyper3[n][[1]][iters], col = 3)
  dev.off()
}

hypercred = t(apply(hyper1, 2, function(x){quantile(x, c(0.025, 0.975))}))

truehyper = c(
  sigma.c = trueParms$sigC,
  d = trueParms$d,
  tau = trueParms$tau,
  theta.phi = trueParms$thePhi,
  theta.alpha = trueParms$theAlp,
  theta.delta = trueParms$theDel,
  sigma.phi = trueParms$sigPhi,
  sigma.alpha = trueParms$sigAlp,
  sigma.delta = trueParms$sigDel,
  pi.alpha = trueParms$piAlp,
  pi.delta = trueParms$piDel
)

truehyper = truehyper[rownames(hypercred)]

hypercred = cbind(hypercred, truehyper)

cover = apply(hypercred, 1, function(x){(x[3] > x[1]) && (x[3] < x[2])})
hypercred = cbind(hypercred, cover)
saveRDS(hypercred, "hypercred.rds")

