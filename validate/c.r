parms1 = read.table("out/parms-chain1.txt", head = T)
parms2 = read.table("out/parms-chain2.txt", head = T)
parms3 = read.table("out/parms-chain3.txt", head = T)

# c

c1 = parms1[, paste("c", 0:5, sep="")]
c2 = parms2[, paste("c", 0:5, sep="")]
c3 = parms3[, paste("c", 0:5, sep="")]

pdf("fig/c-traces-chain1.pdf")
plot(x = 0, y = 0, xlim = c(0, dim(c1)[1]), ylim = c(min(c1), max(c1)), col = "white")
for(p in 1:5)
  lines(x = 1:length(c1[,p]), y = c1[,p], col = p)
dev.off()

pdf("fig/c-traces-chain2.pdf")
plot(x = 0, y = 0, xlim = c(0, dim(c2)[1]), ylim = c(min(c2), max(c2)), col = "white")
for(p in 1:5)
  lines(x = 1:length(c2[,p]), y = c2[,p], col = p)
dev.off()

pdf("fig/c-traces-chain3.pdf")
plot(x = 0, y = 0, xlim = c(0, dim(c3)[1]), ylim = c(min(c3), max(c3)), col = "white")
for(p in 1:5)
  lines(x = 1:length(c3[,p]), y = c3[,p], col = p)
dev.off()