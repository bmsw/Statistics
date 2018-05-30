# NMDS
nmds=metaMDS(matrices,autotransform = T, try=100)
specs=na.omit(nmds$species)
dim(specs)
b=a$FL1.H[as.numeric(rownames(specs))]
b=as.numeric(b)

b[b > 0 & b <= 25] <- 1
b[b > 25 & b <= 50] <- 2
b[b > 50 & b <= 75] <- 3
b
table(b)

plot(nmds, type="n", main="SSC")
points(specs[b==1,], pch=15, cex=.5, col=1)
points(specs[b==2,], pch=15, cex=.5, col=2)
points(specs[b==3,], pch=15, cex=.5, col=3)
points(specs, pch=15, cex=.5, col=b)
text(nmds$points, cex=1.5, labels = c(1:nrow(matrices)), col="blue")


# CCA
nmds=cca(matrices, scale=T )
specs=nmds$CA$v[,1:2]



b=a$FL1.H[as.numeric(rownames(specs))]
b=as.numeric(b)

b[b > 0 & b <= 25] <- 1
b[b > 25 & b <= 50] <- 2
b[b > 50 & b <= 75] <- 3


b[b > 0 & b <= 37] <- 1
b[b > 37 & b <= 75] <- 2

b
table(b)

plot(nmds, type="n", main="SSC")

points(specs[b==1,], pch=15, cex=.5, col=1)
points(specs[b==2,], pch=15, cex=.5, col=2)
points(specs[b==3,], pch=15, cex=.5, col=3)
points(specs, pch=15, cex=.5, col=b)

text(nmds$CA$u[,1:2], cex=1.5, labels = c(1:nrow(matrices)), col="blue")


fit=envfit(nmds, env2[,-1], na.rm=T)
plot(fit, p.max=0.05)