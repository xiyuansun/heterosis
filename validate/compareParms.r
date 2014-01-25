trueParms = readRDS("parms.RDS")
hyper1 = read.table("out/hyper-chain1.txt", head = T)
hyper2 = read.table("out/hyper-chain2.txt", head = T)
hyper3 = read.table("out/hyper-chain3.txt", head = T)

if(!file.exists("fig"))
  dir.create("fig")

for(n in names(hyper)){
  pdf(paste("fig/", n, ".pdf", sep=""))
  plot(x = rep(1:length(hyper1[n][[1]]), 3), y = c(hyper1[n][[1]], hyper2[n][[1]], hyper3[n][[1]]), lty = 1, ylab = n, xlab = "iteration", col = "white")
  lines(x = 1:length(hyper1[n][[1]]), y = hyper1[n][[1]], col = 1)
  lines(x = 1:length(hyper2[n][[1]]), y = hyper2[n][[1]], col = 2)
  lines(x = 1:length(hyper3[n][[1]]), y = hyper3[n][[1]], col = 3)
  dev.off()
}