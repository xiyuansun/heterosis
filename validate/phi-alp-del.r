parms1 = read.table("out/parms-chain1.txt", head = T)
parms2 = read.table("out/parms-chain2.txt", head = T)
parms3 = read.table("out/parms-chain3.txt", head = T)

#phi

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

#alpha

alpha1 = parms1[, paste("alpha", 0:99, sep="")]
alpha2 = parms2[, paste("alpha", 0:99, sep="")]
alpha3 = parms3[, paste("alpha", 0:99, sep="")]

pdf("fig/alpha-traces-chain1.pdf")
plot(x = 0, y = 0, xlim = c(0, dim(alpha1)[1]), ylim = c(min(alpha1), max(alpha1)), col = "white")
for(p in 1:100)
  lines(x = 1:length(alpha1[,p]), y = alpha1[,p], col = p)
dev.off()

pdf("fig/alpha-traces-chain2.pdf")
plot(x = 0, y = 0, xlim = c(0, dim(alpha2)[1]), ylim = c(min(alpha2), max(alpha2)), col = "white")
for(p in 1:100)
  lines(x = 1:length(alpha2[,p]), y = alpha2[,p], col = p)
dev.off()

pdf("fig/alpha-traces-chain3.pdf")
plot(x = 0, y = 0, xlim = c(0, dim(alpha3)[1]), ylim = c(min(alpha3), max(alpha3)), col = "white")
for(p in 1:100)
  lines(x = 1:length(alpha3[,p]), y = alpha3[,p], col = p)
dev.off()

#delta

delta1 = parms1[, paste("delta", 0:99, sep="")]
delta2 = parms2[, paste("delta", 0:99, sep="")]
delta3 = parms3[, paste("delta", 0:99, sep="")]

pdf("fig/delta-traces-chain1.pdf")
plot(x = 0, y = 0, xlim = c(0, dim(delta1)[1]), ylim = c(min(delta1), max(delta1)), col = "white")
for(p in 1:100)
  lines(x = 1:length(delta1[,p]), y = delta1[,p], col = p)
dev.off()

pdf("fig/delta-traces-chain2.pdf")
plot(x = 0, y = 0, xlim = c(0, dim(delta2)[1]), ylim = c(min(delta2), max(delta2)), col = "white")
for(p in 1:100)
  lines(x = 1:length(delta2[,p]), y = delta2[,p], col = p)
dev.off()

pdf("fig/delta-traces-chain3.pdf")
plot(x = 0, y = 0, xlim = c(0, dim(delta3)[1]), ylim = c(min(delta3), max(delta3)), col = "white")
for(p in 1:100)
  lines(x = 1:length(delta3[,p]), y = delta3[,p], col = p)
dev.off()