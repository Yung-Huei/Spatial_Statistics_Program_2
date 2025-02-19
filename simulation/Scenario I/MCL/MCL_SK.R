#simulation
library(expint) #gammainc
library(SpatialTools) #dist1/ dist2/
library(Matrix) #tcrossprod
library(deldir) #deldir
library(stats) #optim/ optimize
library(mvnfast) #dmvn
library(parallel) #makeCluster
library(geoR) #loglik.GRF/ matern
library(foreach) #foreach
library(iterators);library(doParallel) #%dopar%

#spatdiv
library(permute)
library(Rfast)
library(fields)
library(vegan)
colMin = Rfast::colMins

load("simulationMCL.RData")

###############
# Domain grid #
###############
g <- 60
grid1 = grid2 = seq(1/2/g, 1 - 1/2/g, l = g)
G <- expand.grid(grid1, grid2)
grids = matrix(cbind(G$Var1, G$Var2), g^2, 2) 
d <- matrix(,g^2,g^2)
d <- dist1(grids)

###################
# True covariance #
###################
center <- matrix(c(0.25,0.75,0.25,0.75,0.25,0.25,0.75,0.75),4,2)
border = c(0,1,0,1)
K <- dim(center)[1]                                                 #the number of partition
sigmasqvector <- matrix(c(1,1,0.25,0.25),K,1)
alphavector <- matrix(c(16,4,16,4),K,1)
nuvector <- matrix(c(1/2,1/2,1/2,1/2),K,1)
a <- 0.01 
nugget <- 0.1
                
COV <- cov.fun.4(grids, center, border, sigmasqvector, alphavector, nuvector, a)
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
  coords <- matrix(grids[sample400,],N,2)
  Z <- matrix(z[sample400,],N,1)
  
  sigmasqvector_hat <- para[u,2:5]
  alphavector_hat <- para[u,6:9]
  nuvector_hat <- para[u,10:13]
  a_hat <- para[u,14]
  tau_hat <- para[u,15]

  t1 <- proc.time()
  R1 <- cov.fun.4(grids, center, border, sigmasqvector_hat, alphavector_hat, nuvector_hat, a_hat) 
  pred_SK[u,] <- R1[,sample400]%*%solve(R1[sample400,sample400]+tau_hat*diag(N))%*%(Z)
  mspe_SK[u,] <- sum((pred_SK[u,]-y)^2)/g^2
  SPEu_SK[u,] <- (pred_SK[u,]-y)^2
  timeSK[u,] <- (proc.time()-t1)[1:3]
  save.image("MCL_SK.RData")
}
  
round(mean(mspe_SK[,1]),4)
round(sd(mspe_SK[,1])/sqrt(T),4)
round(sapply(1:3, function(j){mean(timeSK[1:T,j])}),4)
round(sapply(1:3, function(j){sd(timeSK[1:T,j])/sqrt(T)}),4)


