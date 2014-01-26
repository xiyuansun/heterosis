trueParms = readRDS("parms.RDS")
hyper1 = read.table("out/hyper-chain1.txt", head = T)
hyper2 = read.table("out/hyper-chain2.txt", head = T)
hyper3 = read.table("out/hyper-chain3.txt", head = T)

if(!file.exists("fig"))
  dir.create("fig")

iters = (length(hyper1[n][[1]]) - 1000):length(hyper1[n][[1]])

for(n in names(hyper1)){
  pdf(paste("fig/", n, ".pdf", sep=""))
  plot(x = rep(iters, 3), y = c(hyper1[n][[1]][iters], hyper2[n][[1]][iters], hyper3[n][[1]][iters]), lty = 1, ylab = n, xlab = "iteration", col = "white")
  lines(x = iters, y = hyper1[n][[1]][iters], col = 1)
  lines(x = iters, y = hyper2[n][[1]][iters], col = 2)
  lines(x = iters, y = hyper3[n][[1]][iters], col = 3)
  dev.off()
}

hypercred = t(apply(hyper1, 2, function(x){quantile(x, c(0.025, 0.975))}))
truehyper = c(
  trueParms$sigC,
  trueParms$d,
  trueParms$tau,
  trueParms$thePhi,
  trueParms$theAlp,
  trueParms$theDel,
  trueParms$sigPhi,
  trueParms$sigAlp,
  trueParms$sigDel,
  trueParms$piAlp,
  trueParms$piDel
)
hypercred = cbind(hypercred, truehyper)

cover = apply(hypercred, 1, function(x){(x[3] > x[1]) && (x[3] < x[2])})
hypercred = cbind(hypercred, cover)
saveRDS(hypercred, "hypercred.rds")

