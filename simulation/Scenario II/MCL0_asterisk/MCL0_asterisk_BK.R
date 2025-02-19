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

load("simulationPAR0.RData")

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
center_T <- matrix(c(0.25,0.75,0.25,0.75,0.25,0.25,0.75,0.75),4,2)
border = c(0,1,0,1)
K <- dim(center_T)[1]                                                 #the number of partition
sigmasqvector <- matrix(c(1,1,0.25,0.25),K,1)
alphavector <- matrix(c(16,4,16,4),K,1)
nuvector <- matrix(c(1/2,1/2,1/2,1/2),K,1)
a <- 0.1 
nugget <- 0.1
                
COV <- cov.fun.4(grids, center_T, border, sigmasqvector, alphavector, nuvector, a)
COV_0.5 <- t(chol(COV))

T <- 100
N <- 1000
pred_BK <- matrix(,T,g^2)
SPEu_BK <- matrix(,T,g^2)
mspe_BK <- matrix(,T,1)
timeBK <- matrix(,T,3)

u <- 1
for( u in 1:T ){
  print(u)
  set.seed(u-1)
  y <- COV_0.5%*%rnorm(g^2,0,1)
  z <- y + sqrt(nugget)*rnorm(g^2,0,1)

  sample400 <- sample(c(1:g^2),N)
  coords <- matrix(grids[sample400,],N,2)
  Z <- matrix(z[sample400,],N,1)
  center <- REP[[u]]

  data_z <- order.fun(as.matrix(Z), as.matrix(coords), center)                  
  grid <- order.fun(as.matrix(coords), as.matrix(coords), center)             
  grid_region <- region.fun(grid, center)                   
  Q <- matrix(,2,K)
  Q[1,] <- cumsum(table(grid_region))-table(grid_region)+1
  Q[2,] <- cumsum(table(grid_region))
  q <- Q[2,]-Q[1,]+1 
  
  sigmasqvector_hat <- para[u,2:5]
  alphavector_hat <- para[u,6:9]
  nuvector_hat <- para[u,10:13]
  a_hat <- para[u,14]
  tau_hat <- para[u,15]

  t1 <- proc.time()   
  Sigma_yii <- lapply(1:K, function (k) {
                 error <- diag(tau_hat, q[k])
                 covariance <- cov.fun.4(grid[Q[1,k]:Q[2,k],], center, border, sigmasqvector_hat, alphavector_hat, nuvector_hat, a_hat)+ error                        
               } )

  Sigma_yii_inverse <- lapply(1:K, function (i) { chol2inv(chol(Sigma_yii[[i]])) })  

  predictsize <- nrow(grids)
  grids_region <- region.fun(grids, center)
  grids_label <- sapply(1:predictsize, function(i){ sum( grids_region[1:i]==grids_region[i] ) })
  pcov <- list()
  pcov <- lapply(1:K, function(i){
            o <- which(grids_region==i)
            if(length(o)==0){
              return(0)
            }else{
              return( cov.fun.3( matrix(grids[o,],length(o),ncol(grids)), 
                                 matrix(grid[Q[1,i]:Q[2,i],],q[i],ncol(grids)), 
                                 center, border, sigmasqvector_hat, alphavector_hat, nuvector_hat, a_hat) )
            }
          })            
 
  pred_BK[u,] <- sapply( 1:predictsize, function (k) { 
                  i <- grids_region[k,]                                   
                  return( pcov[[i]][grids_label[k],] %*% Sigma_yii_inverse[[i]] %*% (data_z[Q[1,i]:Q[2,i]]) )
                 })
  timeBK[u,] <- (proc.time()-t1)[1:3]
  mspe_BK[u,1] <- mean((pred_BK[u,]-y)^2)
  SPEu_BK[u,] <- (pred_BK[u,]-y)^2
  save.image("MCL_BK.RData")
}  
round(mean(mspe_BK[,1]),4)
round(sd(mspe_BK[,1])/sqrt(T),4)
round(sapply(1:3, function(j){mean(timeBK[1:T,j])}),4)
round(sapply(1:3, function(j){sd(timeBK[1:T,j])/sqrt(T)}),4)

  



