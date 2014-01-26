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
  lines(x = iters, y = hyper1[n][[1]], col = 1)
  lines(x = iters, y = hyper2[n][[1]], col = 2)
  lines(x = iters, y = hyper3[n][[1]], col = 3)
  dev.off()
}
