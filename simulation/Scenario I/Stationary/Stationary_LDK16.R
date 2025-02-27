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

load("simulationSTAML.RData")
load("code1e03new.RData")

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

pred_LDK16 <- matrix(,T,g^2)
SPEu_LDK16 <- matrix(,T,g^2)
mspe_LDK16 <- matrix(,T,1)
timeLDK16 <- matrix(,T,3)

u <- 1
for( u in 1:T ){
  print(u)
  set.seed(u-1)
  y <- COV_0.5%*%rnorm(g^2,0,1)
  z <- y + sqrt(nugget)*rnorm(g^2,0,1)

  sample400 <- sample(c(1:g^2),N)
  coords <- matrix(grids[sample400,],N,2)
  Z <- matrix(z[sample400,],N,1)

  ER <- list( center = center,
              border = border,
              Z = Z,
              grid = coords,
              beta0 = 0,
              sigmasq_hat = paraML[u,1],
              alpha_hat = paraML[u,2],
              nu_hat = 0.5,
              tau_hat = paraML[u,3],
              cl = 1000 )

  t1 <- proc.time()
  M <- 16
  sp = seq(1/2/4, 1-1/2/4, l = 4)
  Sp <- expand.grid(sp, sp)
  center_p <- as.matrix(Sp) 
  LDK16 <- STALDK(E=ER, center_p=center_p, grids=grids)

  timeLDK16[u,] <- (proc.time()-t1)[1:3]
  pred_LDK16[u,] <- LDK16$predictor
  mspe_LDK16[u,1] <- mean((pred_LDK16[u,]-y)^2)
  SPEu_LDK16[u,] <- (pred_LDK16[u,]-y)^2
  save.image("STA_LDK16.RData")
}

round(mean(mspe_LDK16[,1]),4)
round(sd(mspe_LDK16[,1])/sqrt(T),4)
round(sapply(1:3, function(j){mean(timeLDK16[1:T,j])}),4)
round(sapply(1:3, function(j){sd(timeLDK16[1:T,j])/sqrt(T)}),4)
