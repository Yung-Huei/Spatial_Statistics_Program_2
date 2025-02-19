#LDK
library(geoR);library(fields);library(MASS)
library(caTools);library(mvnfast);library(psych)
library(ggvoronoi);library(deldir);library(ggplot2);library(SpatialTools);library(cluster)
library(emdbook);library(scoringRules);library(stats);library(Matrix);library(dplyr);library(parallel)
library(GpGp);library(gstat);library(sp)
library(factoextra);library(fpc)
library(nloptr);library(expint);library(stats)

#load("simulationMCL.RData")
load("OURFUN0105.RData")
g <- 60
grid1 = grid2 = seq(1/2/g, 1 - 1/2/g, l = g)
G <- expand.grid(grid1, grid2)
grids = matrix(cbind(G$Var1, G$Var2), g^2, 2) 
d <- matrix(,g^2,g^2)
d <- dist1(grids)

#============================================
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
SN <- 1000

pred_LDK81 <- matrix(,T,g^2)
SPEu_LDK81 <- matrix(,T,g^2)
mspe_LDK81 <- matrix(,T,1)
timeLDK81 <- matrix(,T,3)

sp = seq(1/2/9, 1-1/2/9, l = 9)
Sp <- expand.grid(sp, sp)
center_p <- as.matrix(Sp)  
border = c(0,1,0,1)
vor_dxy_p <- deldir(center_p[,1], center_p[,2], rw = border)              #voronoi tessellation
M <- dim(center_p)[1]
#center_p <- center
#vor_dxy_p <- vor_dxy
#M <- K

u <- 1
for( u in 1:T ){
  print(u)
  set.seed(u-1)
  y <- COV_0.5%*%rnorm(g^2,0,1)
  z <- y + sqrt(nugget)*rnorm(g^2,0,1)

  sample400 <- sample(c(1:g^2),SN)
  coords <- grids[sample400,]
  Z <- z[sample400,]

#  data_z <- order.fun(as.matrix(Z), as.matrix(coords), center)                  
#  grid <- order.fun(as.matrix(coords), as.matrix(coords), center)             
#  grid_region <- region.fun(grid, center)               
#  grid_segment <- d.fun(grid, vor_dxy$dirsgs)    
#  Q <- matrix(,2,K)
#  Q[1,] <- cumsum(table(grid_region))-table(grid_region)+1
#  Q[2,] <- cumsum(table(grid_region))
#  q <- Q[2,]-Q[1,]+1 

data_z_p <- order.fun(as.matrix(Z), as.matrix(coords), center_p)                  
grid_p <- order.fun(as.matrix(coords), as.matrix(coords), center_p) 
grid_p_region <- region.fun(grid_p, center_p)                       
grid_region_p <- region.fun(grid_p, center)
grid_segment_p <- d.fun(grid_p, vor_dxy$dirsgs)
Qmatrix_p <- matrix(,2,M)
Qmatrix_p[1,] <- cumsum(table(grid_p_region))-table(grid_p_region)+1
Qmatrix_p[2,] <- cumsum(table(grid_p_region))
qq <- Qmatrix_p[2,]-Qmatrix_p[1,]+1

#sigmasqvector_hat <- para[u,2:5]
#alphavector_hat <- para[u,6:9]
#nuvector_hat <- para[u,10:13]
#a_hat <- para[u,14]
#tau_hat <- para[u,15]

sigmasqvector_hat <- sigmasqvector
alphavector_hat <- alphavector
nuvector_hat <- nuvector
a_hat <- a
tau_hat <- nugget
alphamatrix_hat <- alpha.fun(as.matrix(alphavector_hat))
numatrix_hat <- nu.fun(as.matrix(nuvector_hat))
rhomatrix_hat <- rho.fun(numatrix_hat, alphamatrix_hat)

t1 <- proc.time()

neighbor_p <- neighbor.fun(center_p, vor_dxy_p)
Neighbor_p <- lapply(1:M, function(i){
  A <- which( neighbor_p[,i]!=0 )
  B <- rbind(as.matrix(vor_dxy_p$dirsgs[A,1:2]),as.matrix(vor_dxy_p$dirsgs[A,3:4]))
  D <- dist1(B)
  del <- c()
  for(j in 1:dim(B)[1]){
    same <- which(D[j,]==0)
    del <-c( del,same[which((same-j)>0)] )
  }
  B <- B[-del,]
  E1 <- dist2(B,as.matrix(vor_dxy_p$dirsgs[,1:2]))
  E2 <- dist2(B,as.matrix(vor_dxy_p$dirsgs[,3:4]))
  F1 <- sapply(1:dim(neighbor_p)[1], function(j){prod(E1[,j])==0})
  F2 <- sapply(1:dim(neighbor_p)[1], function(j){prod(E2[,j])==0})
  relation <- which(F1+F2 != 0)
  G <- sapply(1:M, function(k){sum(neighbor_p[relation,k])})
  return( setdiff(which( G!=0 ),i) )
})


Sigma_yii <- lapply(1:M, function (k) {
                                       error <- diag(tau_hat, qq[k])
                                       covariance <- cov.fun(grid_p[Qmatrix_p[1,k]:Qmatrix_p[2,k],], 
neighbormatrix, 
as.matrix(grid_segment_p[Qmatrix_p[1,k]:Qmatrix_p[2,k],]), 
as.matrix(grid_region_p[Qmatrix_p[1,k]:Qmatrix_p[2,k],]), 
sigmasqvector_hat, rhomatrix_hat, alphamatrix_hat, numatrix_hat, a_hat)+ error                        
                                     } )

Sigma_yij <- list()
for(j in 1:(M-1)){
  N <- length(Neighbor_p[[j]])
  A1=A <- Neighbor_p[[j]]
  for(i in 1:N){
    A <- unique(c(A,Neighbor_p[[A1[i]]])) 
  }
  A <- sort(A)   
  Sigma_yij[[j]] <- lapply((j+1):M , function (i) {
    if(sum(i == A)){ 
      covariance <- cov.fun.2(qq[j],grid_p[c(Qmatrix_p[1,j]:Qmatrix_p[2,j],Qmatrix_p[1,i]:Qmatrix_p[2,i]),], neighbormatrix, 
as.matrix(grid_segment_p[c(Qmatrix_p[1,j]:Qmatrix_p[2,j],Qmatrix_p[1,i]:Qmatrix_p[2,i]),]), 
as.matrix(grid_region_p[c(Qmatrix_p[1,j]:Qmatrix_p[2,j],Qmatrix_p[1,i]:Qmatrix_p[2,i]),]), 
sigmasqvector_hat, rhomatrix_hat, alphamatrix_hat, numatrix_hat, a_hat)
    }else{
      return(0)
    }                      
  })
}

Sigma_yii_inverse <- lapply(1:M, function (i) {chol2inv(chol(Sigma_yii[[i]]))})  

predictsize <- nrow(grids)
predgrids_region <- region.fun(grids, center_p)
pred_LDK81[u,] <- sapply( 1:predictsize, function (k) {
   u <- predgrids_region[k,]
   N <- length(Neighbor_p[[u]])+1
   A <- sort(c(u,Neighbor_p[[u]]))
   pcov <- lapply(1:N, function(i){cov.fun.2(1,
rbind(matrix(grids[k,],1,2),grid_p[Qmatrix_p[1,A[i]]:Qmatrix_p[2,A[i]],]), 
neighbormatrix, 
as.matrix(rbind(grids_segment[k,],grid_segment_p[Qmatrix_p[1,A[i]]:Qmatrix_p[2,A[i]],])),
rbind(as.matrix(grids_region[k,]),as.matrix(grid_region_p[Qmatrix_p[1,A[i]]:Qmatrix_p[2,A[i]],])), 
sigmasqvector_hat, rhomatrix_hat, alphamatrix_hat, numatrix_hat, a_hat) })                                    
   bi <- lapply(1:N, function (i) { crossprod(t(pcov[[i]]) , Sigma_yii_inverse[[A[i]]]) } )
   hat_z_local <- sapply(1:N, function (i) { crossprod(t(bi[[i]]) , data_z_p[Qmatrix_p[1,A[i]]:Qmatrix_p[2,A[i]]]) } )
   di <- as.matrix(sapply(1:N, function (i) { tcrossprod(bi[[i]] , pcov[[i]]) } ))
   MM <- matrix(, N, N)
   for (i in 1:N) {
     for (j in 1:N) {
       if(i==j){ 
          MM[i,j] <- crossprod( t(bi[[i]]) , tcrossprod(Sigma_yii[[A[i]]] , bi[[j]]) )
       }else if (A[i]>A[j]){
          MM[i,j] <- crossprod( t(bi[[i]]) , tcrossprod(t(Sigma_yij[[A[j]]][[A[i]-A[j]]]) , bi[[j]]) )
       }else if (A[i]<A[j]){
          MM[i,j] <- crossprod( t(bi[[i]]) , tcrossprod(Sigma_yij[[A[i]]][[A[j]-A[i]]] , bi[[j]]) )
       }  
     }
   }
   return( crossprod(di, chol2inv(chol(MM))) %*% hat_z_local )
})

timeLDK81[u,] <- (proc.time()-t1)[1:3]
mspe_LDK81[u,1] <- mean((pred_LDK81[u,]-y)^2)
SPEu_LDK81[u,] <- (pred_LDK81[u,]-y)^2
save.image("TRUE_LDK81.RData")
}

round(mean(mspe_LDK81[,1]),4)
round(sd(mspe_LDK81[,1])/sqrt(T),4)
round(mean(timeLDK81[,1]),4)
round(sd(timeLDK81[,1])/sqrt(T),4)
write.csv(mspe_LDK81, "LDK81.csv", row.names=FALSE)
write.csv(timeLDK81, "timeLDK81.csv", row.names=FALSE)


