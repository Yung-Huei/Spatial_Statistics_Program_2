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
library(psych) #tr

#spatdiv
library(permute)
library(Rfast)
library(fields)
library(vegan)
colMin = Rfast::colMins

load("code1e03new.RData")
#load("spatdiv.RData")

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
sigmasqvector <- matrix(c(1,1,1,1),K,1)
alphavector <- matrix(c(4,4,4,4),K,1)
nuvector <- matrix(c(1/2,1/2,1/2,1/2),K,1)
a <- 0 
nugget <- 0.1
                
COV <- cov.fun.4(grids, center, border, sigmasqvector, alphavector, nuvector, a)
COV_0.5 <- t(chol(COV))

T <- 100
N <- 1000
para <- matrix(,T,3*K+4) 
steptime <- matrix(,T,4)
F_loss <- matrix(,T,1)
KL_loss <- matrix(,T,1) 

u <- 1
for( u in 1:T ){
  print(u)
  set.seed(u-1)
  y <- COV_0.5%*%rnorm(g^2,0,1)
  z <- y + sqrt(nugget)*rnorm(g^2,0,1)

  sample400 <- sample(c(1:g^2),N)
  coords <- matrix(grids[sample400,],N,2)
  Z <- matrix(z[sample400,],N,1)

  ER <- MCL_exp(grid=coords, Z=Z, center=center, border=border)

  para[u,] <- c(0, ER$sigmasqvector_hat, ER$alphavector_hat, ER$nuvector_hat, ER$a_hat, ER$tau_hat, ER$cl)
  steptime[u,] <- ER$time[,1]
  B <- cov.fun.4(grids, center, border, ER$sigmasqvector_hat, ER$alphavector_hat, ER$nuvector_hat, ER$a_hat) 
  C <- COV
  A <- B - C
  F_loss[u,] <- sqrt(tr(A%*%t(A)))
  KL_loss[u,] <- 0.5*(tr(solve(B)%*%C)-3600+determinant(B)$modulus-determinant(C)$modulus)
  
  save.image("simulationMCL.RData")
  if(u==10| u==50){
  print(round(sapply(1:(3*K+4), function(j){mean(para[1:u,j])}),4))
  print(round(sapply(1:(3*K+4), function(j){sd(para[1:u,j])/sqrt(u)}),4))
  print(round(sapply(1:4, function(j){mean(steptime[1:u,j])}),4))
  print(round(sapply(1:4, function(j){sd(steptime[1:u,j])/sqrt(u)}),4))
  print(round(mean(F_loss[1:u,1]),4))
  print(round(sd(F_loss[1:u,1])/sqrt(u),4))
  print(round(mean(KL_loss[1:u,1]),4))
  print(round(sd(KL_loss[1:u,1])/sqrt(u),4))
  }
}

round(sapply(1:(3*K+4), function(j){mean(para[1:T,j])}),4)
round(sapply(1:(3*K+4), function(j){sd(para[1:T,j])}),4)
round(sapply(1:4, function(j){mean(steptime[1:T,j])}),4)
round(sapply(1:4, function(j){sd(steptime[1:T,j])/sqrt(T)}),4)
round(mean(F_loss[,1]),4)
round(sd(F_loss[,1]),4)
round(mean(KL_loss[,1]),4)
round(sd(KL_loss[,1]),4)
