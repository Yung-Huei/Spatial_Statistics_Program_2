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
sigmasqvector <- matrix(c(1,1,0.25,0.25),K,1)
alphavector <- matrix(c(16,4,16,4),K,1)
nuvector <- matrix(c(1/2,1/2,1/2,1/2),K,1)
a <- 0.1 
nugget <- 0.1
                
COV <- cov.fun.4(grids, center, border, sigmasqvector, alphavector, nuvector, a)
COV_0.5 <- t(chol(COV))

T <- 100
N <- 1000
paraML <- matrix(,T,3) 
timeML <- matrix(,T,3)
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
  
  t1 <- proc.time()
  like <- function(par){
           loglik.GRF2(coords = coords, data = Z, 
           obj.model = NULL, cov.model = "exp", cov.pars=c(par[1],1/par[2]), nugget = par[3], 
           kappa = 0.5, lambda = 1, psiR = 1, psiA = 0, trend = "cte", 
           method.lik = "ML", compute.dists = TRUE, realisations = NULL) }
  E <- optim(par = c(mean(sigmasqvector), mean(alphavector), nugget), like, lower = c(rep(1e-06,3)), upper = c(rep(Inf,3)), 
           method = "L-BFGS-B", control = list(fnscale=-1))
  timeML[u,] <- (proc.time()-t1)[1:3]
  paraML[u,] <- E$par  

  B <- E$par[1]*matern(d, 1/E$par[2], 0.5) 
  C <- COV 
  A <- B - C
  F_loss[u,] <- sqrt(tr(A%*%t(A)))
  KL_loss[u,] <- 0.5*(tr(solve(B)%*%C)-3600+determinant(B)$modulus-determinant(C)$modulus)

  save.image("simulationSTAML.RData")
  if(u==10| u==50){
  print(round(sapply(1:3, function(j){mean(paraML[1:u,j])}),4))
  print(round(sapply(1:3, function(j){sd(paraML[1:u,j])/sqrt(u)}),4))
  print(round(mean(F_loss[1:u,1]),4))
  print(round(sd(F_loss[1:u,1])/sqrt(u),4))
  print(round(mean(KL_loss[1:u,1]),4))
  print(round(sd(KL_loss[1:u,1])/sqrt(u),4))
  print(round(sapply(1:3, function(j){mean(timeML[1:u,j])}),4))
  print(round(sapply(1:3, function(j){sd(timeML[1:u,j])/sqrt(u)}),4))
  }
}
  
round(sapply(1:3, function(j){mean(paraML[1:T,j])}),4)
round(sapply(1:3, function(j){sd(paraML[1:T,j])/sqrt(u)}),4)
round(mean(F_loss[1:T,1]),4)
round(sd(F_loss[1:T,1])/sqrt(u),4)
round(mean(KL_loss[1:T,1]),4)
round(sd(KL_loss[1:T,1])/sqrt(u),4)
round(sapply(1:3, function(j){mean(timeML[1:T,j])}),4)
round(sapply(1:3, function(j){sd(timeML[1:T,j])/sqrt(u)}),4)



