#LDK
library(geoR);library(fields);library(MASS)
library(caTools);library(mvnfast);library(psych)
library(ggvoronoi);library(deldir);library(ggplot2);library(SpatialTools);library(cluster)
library(emdbook);library(scoringRules);library(stats);library(Matrix);library(dplyr);library(parallel)
library(GpGp);library(gstat);library(sp)
library(factoextra);library(fpc)
library(nloptr);library(expint);library(stats)


load("OURFUN0105.RData")
g <- 60
grid1 = grid2 = seq(1/2/g, 1 - 1/2/g, l = g)
G <- expand.grid(grid1, grid2)
grids = matrix(cbind(G$Var1, G$Var2), g^2, 2) 
d <- matrix(,g^2,g^2)
d <- dist1(grids)

#parameter settings ============================================
center <- matrix(c(0.25,0.75,0.25,0.75,0.25,0.25,0.75,0.75),4,2)
border = c(0,1,0,1)
vor_dxy <- deldir(center[,1], center[,2], rw = border)              #voronoi tessellation
K <- dim(center)[1]                                                 #the number of partition
neighbormatrix <- neighbor.fun(center, vor_dxy)
grids_segment <- d.fun(grids, vor_dxy$dirsgs) 
grids_region <- region.fun(grids, center) 
sigmasqvector <- matrix(c(1,1,0.25,0.25),K,1)
alphavector <- matrix(c(16,4,16,4),K,1)
nuvector <- matrix(c(1/2,1/2,1/2,1/2),K,1)
nugget <- 0.1
alphamatrix <- alpha.fun(alphavector)
numatrix <- nu.fun(nuvector)
rhomatrix <- rho.fun(numatrix, alphamatrix)
a <- 0.01                 
COV <- cov.fun(grids, neighbormatrix, grids_segment, grids_region, sigmasqvector, rhomatrix, alphamatrix, numatrix, a)
COV_0.5 <- t(chol(COV))
T <- 100
N <- 1000

pred_SK <- matrix(,T,g^2)
SPEu_SK <- matrix(,T,g^2)
mspe_SK <- matrix(,T,1)
timeSK <- matrix(,T,3)

u <- 1
for( u in 1:T ){
  print(u)
  set.seed(u-1)
  y <- COV_0.5%*%rnorm(g^2,0,1)
  z <- y + sqrt(nugget)*rnorm(g^2,0,1)

  sample400 <- sample(c(1:g^2),N)
  coords <- grids[sample400,]
  Z <- z[sample400,]

  data_z <- order.fun(as.matrix(Z), as.matrix(coords), center)                  
  grid <- order.fun(as.matrix(coords), as.matrix(coords), center)             
  grid_region <- region.fun(grid, center)               
  grid_segment <- d.fun(grid, vor_dxy$dirsgs)    
  Q <- matrix(,2,K)
  Q[1,] <- cumsum(table(grid_region))-table(grid_region)+1
  Q[2,] <- cumsum(table(grid_region))
  q <- Q[2,]-Q[1,]+1 
  
sigmasqvector_hat <- sigmasqvector
alphavector_hat <- alphavector
nuvector_hat <- nuvector
a_hat <- a
tau_hat <- nugget
alphamatrix_hat <- alpha.fun(as.matrix(alphavector_hat))
numatrix_hat <- nu.fun(as.matrix(nuvector_hat))
rhomatrix_hat <- rho.fun(numatrix_hat, alphamatrix_hat)

t1 <- proc.time()
R1 <- cov.fun(grids, neighbormatrix, grids_segment, grids_region, sigmasqvector_hat, rhomatrix_hat, alphamatrix_hat, numatrix_hat, a_hat) 
pred_SK[u,] <- R1[,sample400]%*%solve(R1[sample400,sample400]+tau_hat*diag(N))%*%(Z)
mspe_SK[u,] <- sum((pred_SK[u,]-y)^2)/g^2
#SPEu_SK[u,] <- (pred_SK[u,]-y)^2
timeSK[u,] <- (proc.time()-t1)[1:3]
save.image("TRUE_SK.RData")
}
  
round(mean(mspe_SK[,1]),4)
round(sd(mspe_SK[,1])/sqrt(T),4)
round(sapply(1:3, function(j){mean(timeSK[1:T,j])}),4)
round(sapply(1:3, function(j){sd(timeSK[1:T,j])/sqrt(T)}),4)
write.csv(mspe_SK, "SK.csv", row.names=FALSE)
write.csv(timeSK, "timeSK.csv", row.names=FALSE)


