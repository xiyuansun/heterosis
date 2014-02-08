parms1 = read.table("out/parms-chain1.txt", head = T)
parms2 = read.table("out/parms-chain2.txt", head = T)
parms3 = read.table("out/parms-chain3.txt", head = T)

phi1 = parms1[, paste("phi", 0:99, sep="")]
phi2 = parms2[, paste("phi", 0:99, sep="")]
phi3 = parms3[, paste("phi", 0:99, sep="")]

pdf("fig/phi-traces-chain1.pdf")
plot(x = 0, y = 0, xlim = c(0, dim(phi1)[1]), ylim = c(min(phi1), max(phi1)), col = "white")
for(p in 1:100)
  lines(x = 1:length(phi1[,p]), y = phi1[,p], col = p)
dev.off()

pdf("fig/phi-traces-chain2.pdf")
plot(x = 0, y = 0, xlim = c(0, dim(phi2)[1]), ylim = c(min(phi2), max(phi2)), col = "white")
for(p in 1:100)
  lines(x = 1:length(phi2[,p]), y = phi2[,p], col = p)
dev.off()

pdf("fig/phi-traces-chain3.pdf")
plot(x = 0, y = 0, xlim = c(0, dim(phi3)[1]), ylim = c(min(phi3), max(phi3)), col = "white")
for(p in 1:100)
  lines(x = 1:length(phi3[,p]), y = phi3[,p], col = p)
dev.off()