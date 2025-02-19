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
sigmasqvector <- matrix(c(1,1,1,1),K,1)
alphavector <- matrix(c(4,4,4,4),K,1)
nuvector <- matrix(c(1/2,1/2,1/2,1/2),K,1)
a <- 0 
nugget <- 0.1
                
COV <- cov.fun.4(grids, center, border, sigmasqvector, alphavector, nuvector, a)
COV_0.5 <- t(chol(COV))

T <- 100
N <- 1000
pred_DK4 <- matrix(,T,g^2)
SPEu_DK4 <- matrix(,T,g^2)
mspe_DK4 <- matrix(,T,1)
timeDK4 <- matrix(,T,3)

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
              sigmasqvector_hat = matrix(para[u,2:5],K,1),
              alphavector_hat = matrix(para[u,6:9],K,1),
              nuvector_hat = matrix(para[u,10:13],K,1),
              a_hat = para[u,14],
              tau_hat = para[u,15],
              cl = 1000 )

  t1 <- proc.time()
  M <- K
  DK4 <- DK(E=ER, center_p=ER$center, grids=grids)

  timeDK4[u,] <- (proc.time()-t1)[1:3]
  pred_DK4[u,] <- DK4$predictor
  mspe_DK4[u,1] <- mean((pred_DK4[u,]-y)^2)
  SPEu_DK4[u,] <- (pred_DK4[u,]-y)^2
  save.image("MCL_DK4.RData")
}

round(mean(mspe_DK4[,1]),4)
round(sd(mspe_DK4[,1])/sqrt(T),4)
round(sapply(1:3, function(j){mean(timeDK4[1:T,j])}),4)
round(sapply(1:3, function(j){sd(timeDK4[1:T,j])/sqrt(T)}),4)
