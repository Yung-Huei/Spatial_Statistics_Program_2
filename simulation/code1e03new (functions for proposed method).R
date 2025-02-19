#Date: 20240620
#Title: Fast Spatial Prediction for Nonstationary Processes 
#       with a Divide-and-Conquer Strategy
#Main Function: 
#Input:
#Output:   
##############################################################

###### load packages ######

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
library(flexclust) #dist2
#library(dplyr)
#library(spaMM)
#library(magrittr)

###### subfunctions ######

nu.fun <- function (nuvector){

  partition <- nrow(nuvector)
  numatrix <- matrix(,partition, partition) 
  for (j in 1:partition) {
      numatrix[,j] <- sapply(1:partition, function(i){( nuvector[i] + nuvector[j] ) / 2 } )
  }
  return(numatrix)
}


alpha.fun <- function (alphavector){

  partition <- nrow(alphavector)
  alphamatrix <- matrix(,partition, partition) 
  for (j in 1:partition) {
      alphamatrix[,j] <- sapply(1:partition, function(i){sqrt( ( alphavector[i]^2 + alphavector[j]^2 ) / 2 ) } )
  }
  return(alphamatrix)
}


rho.fun <- function (numatrix, alphamatrix){

  partition <- nrow(numatrix)
  rhomatrix <- matrix(,partition, partition)    
  for (k in 1:partition) {
      rhomatrix[,k] <- sapply(1:partition, function(i){
        alphamatrix[i,i]^numatrix[i,i] * alphamatrix[k,k]^numatrix[k,k] / ( alphamatrix[i,k]^(2*numatrix[i,k]) ) * 
        gammainc(numatrix[i,k], 0) / ( sqrt (gammainc(numatrix[i,i], 0) * gammainc(numatrix[k,k], 0)) )  } )
  }
  return(rhomatrix)
}


region.fun <- function (grid, center){

  gridsize <- nrow(grid)
  grid_center <- dist2(grid, center)
  region <- as.matrix(sapply(1:gridsize, function (i) {order(grid_center[i,])[1]} ))
  return(region)
}


order.fun <- function (Y, grid, center){

  Y_region <- region.fun(grid, center)
  Y_order <- Y[order(Y_region),]
  return(as.matrix(Y_order))
}


#Output:neighbor(分割線段*區塊數)
#               (第v個分割線段是否是第u塊邊界，"是"則非0且等於某常數c，代表分割第u塊和第c塊)
neighbor.fun <- function (center, border){

  partition <- nrow(center)
  voronoi <- deldir(center[,1], center[,2], rw = border)$dirsgs    
  neighbor <- matrix(0, nrow(voronoi), partition)
  for (u in 1:partition) {
      for (v in 1:dim(voronoi)[1]) {
          if ((voronoi$ind1[v] - u) * (voronoi$ind2[v] - u) == 0) {
             neighbor[v,u] <- voronoi$ind1[v] + voronoi$ind2[v] - u
          }
      }
  }
  return(neighbor)
}


d.fun <- function (grid, center, border){
  
  dirsgs <- deldir(center[,1], center[,2], rw = border)$dirsgs  
  grid_dirsgs <- matrix(,nrow(grid), nrow(dirsgs))
  for (j in 1:nrow(dirsgs)) {
      for (i in 1:dim(grid)[1]) { 
          A <- c(dirsgs[j,1], dirsgs[j,2])
          B <- c(dirsgs[j,3], dirsgs[j,4])
          P <- grid[i,]         
          AB <- B - A
          AP <- P - A
          PB <- B - P
          r <- as.numeric( ifelse(AB%*%AB == 0, 0, AB %*% AP / AB %*% AB) )
          if (0 < r && r < 1) {
             C <- r * AB + A
             PC <- C - P
             grid_dirsgs[i,j] <- sqrt( PC %*% PC )
          }else if (r <= 0) {
                   grid_dirsgs[i,j] <- sqrt( AP %*% AP )
          }else{
               grid_dirsgs[i,j] <- sqrt( PB %*% PB )
          }
      }
  }
  return(grid_dirsgs)
}


weight.fun1 <- function (grid, center, border, a_parameter){
 
  neighbormatrix <- neighbor.fun(center, border)
  N <- nrow(grid)
  K <- nrow(center)
  weight <- matrix(, N, K)  #權重函數 (格點數量*區塊數)
  for (u in 1:N) {
      grid_segment <- d.fun(matrix(grid[u,],1,ncol(grid)), center, border)
      grid_region <- region.fun(matrix(grid[u,],1,ncol(grid)), center)
      for (v in 1:K) {
          if(a_parameter == 0){
             weight[u,v] <- 1*(grid_region == v)
          }else if(grid_region == v) {
             weight[u,v] <- 1
          }else{
             weight[u,v] <- 1 - min( ( min( grid_segment[which( neighbormatrix[,v] != 0 )] ) / a_parameter ), 1 ) 
          }         
      }
  }
  return(weight)
}


weight.fun2 <- function (grid, center, border, rhomatrix, a_par){

  N <- nrow(grid)
  K <- nrow(center)
  weight <- matrix(,N, K)
  w <- weight.fun1(grid, center, border, a_par)
  for(u in 1:N){   
     A <- 0
     for(i in 1:K){
       for(j in 1:K){ 
         A <- A + rhomatrix[i,j] * w[u,i] * w[u,j]
       }
     }
     for(v in 1:K){ 
         weight[u,v] <- w[u,v] / sqrt(A)
     }
  }
  return(weight)
}


#old cov.fun
cov.fun <- function (grid, center, border, sigmasqvector, alphavector, nuvector, a_par){

  N <- nrow(grid)
  K <- nrow(center) 
  alphamatrix <- alpha.fun(matrix(alphavector,K,1)) 
  numatrix <- nu.fun(matrix(nuvector,K,1))
  rhomatrix <- rho.fun(numatrix, alphamatrix)

  grid_grid <- as.matrix(dist1(grid))  
  weight <- weight.fun2(grid, center, border, rhomatrix, a_par)
  cov <- matrix(0, N, N)
  NK <- expand.grid(1:K, 1:K)
  results <- lapply(1:K^2, function(i){
               k <- NK[i,1]
               l <- NK[i,2]
               M <- matern(grid_grid, phi = 1/alphamatrix[k,l], kappa = numatrix[k,l])
               return( rhomatrix[k,l] * sqrt( sigmasqvector[k] ) * sqrt( sigmasqvector[l] ) * tcrossprod(weight[,k],weight[,l]) * M )
             })                             
  cov <- Reduce('+',results)
  return(cov)
}


#old cov.fun
cov.fun.2 <- function (cut, grid, center, border, sigmasqvector, alphavector, nuvector, a_par){

  N <- nrow(grid)
  K <- nrow(center) 
  alphamatrix <- alpha.fun(matrix(alphavector,K,1)) 
  numatrix <- nu.fun(matrix(nuvector,K,1))
  rhomatrix <- rho.fun(numatrix, alphamatrix)

  grid_grid <- as.matrix( dist2(matrix(grid[1:cut,], cut, 2), matrix(grid[-(1:cut),], N-cut, 2)) )  
  weight <- weight.fun2(grid, center, border, rhomatrix, a_par)
  cov <- matrix(0, cut, N-cut)
  NK <- expand.grid(1:K, 1:K)
  results <- lapply(1:K^2, function(i){
               k <- NK[i,1]
               l <- NK[i,2]  
               M <- matern(grid_grid, phi = 1/alphamatrix[k,l], kappa = numatrix[k,l])
               rhomatrix[k,l] * sqrt( sigmasqvector[k] ) * sqrt( sigmasqvector[l] ) * 
                 matrix(weight[1:cut,k],cut,1) %*% matrix(weight[-(1:cut),l],1,N-cut) * M                            
             })                             
  cov <- Reduce('+',results)
  return(cov)
}


cov.fun.3 <- function (grid1, grid2, center, border, sigmasqvector, alphavector, nuvector, a_par) {
  
  N1 <- nrow(grid1)
  N2 <- nrow(grid2)
  K <- nrow(center)  
  alphamatrix <- alpha.fun(matrix(alphavector, K, 1))
  numatrix <- nu.fun(matrix(nuvector, K, 1))
  rhomatrix <- rho.fun(numatrix, alphamatrix)
  
  neighbor <- neighbor.fun(center, border)
  vor_dxy <- deldir(center[, 1], center[, 2], rw = border)$dirsgs  
  Neighbor <- lapply(1:K, function(i) {
    A <- which(neighbor[, i] != 0)
    B <- rbind(as.matrix(vor_dxy[A, 1:2]), as.matrix(vor_dxy[A, 3:4]))
    D <- dist1(B)
    del <- c()
    for (j in 1:dim(B)[1]) {
      same <- which(D[j, ] == 0)
      del <- c(del, same[which((same - j) > 0)])
    }
    if(length(del)!=0){
        B <- B[-del, ]
    }
    E1 <- dist2(B, as.matrix(vor_dxy[, 1:2]))
    E2 <- dist2(B, as.matrix(vor_dxy[, 3:4]))
    F1 <- sapply(1:dim(neighbor)[1], function(j) { prod(E1[, j]) == 0 })
    F2 <- sapply(1:dim(neighbor)[1], function(j) { prod(E2[, j]) == 0 })
    relation <- which(F1 + F2 != 0)
    G <- sapply(1:K, function(k) { sum(neighbor[relation, k]) })
    return(setdiff(which(G != 0), i))
  })
  
  grid1_grid2 <- matrix(dist2(grid1, grid2), N1, N2) 
  weight1 <- weight.fun2(grid1, center, border, rhomatrix, a_par)
  weight2 <- weight.fun2(grid2, center, border, rhomatrix, a_par)
  grid1_region <- region.fun(grid1, center)
  grid2_region <- region.fun(grid2, center)
  NK <- expand.grid(1:K, 1:K)
  cov <- matrix(0, N1, N2)
  for(i in 1:K^2){
      m <- NK[i,1]
      n <- NK[i,2]
      o1 <- which(grid1_region==m)
      o2 <- which(grid2_region==n) 
      neighbors_m <- c(m, Neighbor[[m]])
      neighbors_n <- c(n, Neighbor[[n]])
      neighbors_mn <- expand.grid(neighbors_m, neighbors_n)
      results <- lapply(1:nrow(neighbors_mn), function(i){
                   k <- neighbors_mn[i,1]
                   l <- neighbors_mn[i,2]
                   M <- matern(matrix(grid1_grid2[o1,o2],length(o1),length(o2)), 
                               phi = 1/alphamatrix[k,l], kappa = numatrix[k,l])              
                   rhomatrix[k,l] * sqrt( sigmasqvector[k] ) * sqrt( sigmasqvector[l] ) * 
                   matrix(weight1[o1,k],length(o1),1) %*% matrix(weight2[o2,l],1,length(o2)) * M })
      temp <- Reduce('+',results)
      cov[o1,o2] <- temp                            
  }                             
  return(cov)
}


cov.fun.4 <- function (grid1, center, border, sigmasqvector, alphavector, nuvector, a_par) {
  
  N1 <- nrow(grid1)
  K <- nrow(center)  
  alphamatrix <- alpha.fun(matrix(alphavector, K, 1))
  numatrix <- nu.fun(matrix(nuvector, K, 1))
  rhomatrix <- rho.fun(numatrix, alphamatrix)
  
  neighbor <- neighbor.fun(center, border)
  vor_dxy <- deldir(center[, 1], center[, 2], rw = border)$dirsgs  
  Neighbor <- lapply(1:K, function(i) {
    A <- which(neighbor[, i] != 0)
    B <- rbind(as.matrix(vor_dxy[A, 1:2]), as.matrix(vor_dxy[A, 3:4]))
    D <- dist1(B)
    del <- c()
    for (j in 1:dim(B)[1]) {
      same <- which(D[j, ] == 0)
      del <- c(del, same[which((same - j) > 0)])
    }
    if(length(del)!=0){
        B <- B[-del, ]
    }
    E1 <- dist2(B, as.matrix(vor_dxy[, 1:2]))
    E2 <- dist2(B, as.matrix(vor_dxy[, 3:4]))
    F1 <- sapply(1:dim(neighbor)[1], function(j) { prod(E1[, j]) == 0 })
    F2 <- sapply(1:dim(neighbor)[1], function(j) { prod(E2[, j]) == 0 })
    relation <- which(F1 + F2 != 0)
    G <- sapply(1:K, function(k) { sum(neighbor[relation, k]) })
    return(setdiff(which(G != 0), i))
  })
  
  grid1_grid2 <- matrix(dist1(grid1), N1, N1)  
  weight1 <- weight.fun2(grid1, center, border, rhomatrix, a_par)
  grid1_region <- region.fun(grid1, center)
  NK <- expand.grid(1:K, 1:K)
  cov <- matrix(0, N1, N1)
  for(i in 1:K^2){
      m <- NK[i,1]
      n <- NK[i,2]
      o1 <- which(grid1_region==m)
      o2 <- which(grid1_region==n) 
      neighbors_m <- c(m, Neighbor[[m]])
      neighbors_n <- c(n, Neighbor[[n]])
      neighbors_mn <- expand.grid(neighbors_m, neighbors_n)
      results <- lapply(1:nrow(neighbors_mn), function(i){
                   k <- neighbors_mn[i,1]
                   l <- neighbors_mn[i,2]
                   M <- matern(matrix(grid1_grid2[o1,o2],length(o1),length(o2)), 
                               phi = 1/alphamatrix[k,l], kappa = numatrix[k,l])              
                   rhomatrix[k,l] * sqrt( sigmasqvector[k] ) * sqrt( sigmasqvector[l] ) * 
                   matrix(weight1[o1,k],length(o1),1) %*% matrix(weight1[o2,l],1,length(o2)) * M })
      temp <- Reduce('+',results)
      cov[o1,o2] <- temp                            
  }                             
  return(cov)
}


#beta0=0
loglik.GRF2 <- function (geodata, coords = geodata$coords, data = geodata$data, 
    obj.model = NULL, cov.model = "exp", cov.pars, nugget = 0, 
    kappa = 0.5, lambda = 1, psiR = 1, psiA = 0, trend = "cte", 
    method.lik = "ML", compute.dists = TRUE, realisations = NULL){

    if (!is.null(obj.model)) {
        if (!is.null(obj.model$cov.model)) 
            cov.model <- obj.model$cov.model
        if (!is.null(obj.model$cov.pars)) 
            cov.pars <- obj.model$cov.pars
        if (!is.null(obj.model$nugget)) 
            nugget <- obj.model$nugget
        if (!is.null(obj.model$kappa)) 
            kappa <- obj.model$kappa
        if (!is.null(obj.model$lambda)) 
            lambda <- obj.model$lambda
        if (!is.null(obj.model$psiR)) 
            psiR <- obj.model$psiR
        if (!is.null(obj.model$psiA)) 
            psiA <- obj.model$psiA
        if (!is.null(obj.model$trend)) 
            trend <- eval(obj.model$trend)
    }
    sigmasq <- cov.pars[1]
    phi <- cov.pars[2]
    if (method.lik == "REML" | method.lik == "reml" | method.lik == 
        "rml") 
        method.lik <- "RML"
    if (method.lik == "ML" | method.lik == "ml") 
        method.lik <- "ML"
    if (is.null(realisations)) 
        realisations <- as.factor(rep(1, length(data)))
    else realisations <- as.factor(realisations)
    nrep <- length(levels(realisations))
    if (kappa < 1e-04) 
        return(-(.Machine$double.xmax^0.5))
    if ((nugget + sigmasq) < 1e-16) 
        return(-(.Machine$double.xmax^0.5))
    if (missing(geodata)) 
        xmat <- unclass(trend.spatial(trend = trend, geodata = list(coords = coords, 
            data = data)))
    else xmat <- unclass(trend.spatial(trend = trend, geodata = geodata))
    if (nrow(xmat) != nrow(coords)) 
        stop("coords and trend have incompatible sizes")
    beta.size <- ncol(xmat)
    xmat <- split(as.data.frame(unclass(xmat)), realisations)
    xmat <- lapply(xmat, as.matrix)
    vecdist <- function(x) {
        as.vector(dist(x))
    }
    if (psiR != 1 | psiA != 0) {
        coords.c <- coords.aniso(coords, aniso.pars = c(psiA, 
            psiR))
        .likGRF.dists.vec <- lapply(split(as.data.frame(coords.c), 
            as.factor(realisations)), vecdist)
    }
    else if (compute.dists) 
        .likGRF.dists.vec <- lapply(split(as.data.frame(coords), 
            as.factor(realisations)), vecdist)
    z <- data
    if (abs(lambda - 1) < 1e-04) 
        log.jacobian <- 0
    else {
        if (any(z <= 0)) 
            stop("Transformation not allowed for zero or negative data")
        data <- z^(lambda - 1)
        if (any(data <= 0)) 
            log.jacobian <- log(prod(data))
        else log.jacobian <- sum(log(data))
        data <- NULL
        if (abs(lambda) < 1e-04) 
            data <- log(z)
        else data <- ((z^lambda) - 1)/lambda
    }
    data <- split(data, as.factor(realisations))
    sumnegloglik <- 0
    for (i in 1:nrep) {
        n <- length(data[[1]])
        if ((phi < 1e-16) | (sigmasq < 1e-16)) {
            V <- list(varcov = diag(x = (nugget + sigmasq), n), 
                log.det.to.half = (n/2) * log(nugget + sigmasq))
        }
        else {
            V <- varcov.spatial(dists.lowertri = .likGRF.dists.vec[[i]], 
                cov.model = cov.model, kappa = kappa, nugget = nugget, 
                cov.pars = c(sigmasq, phi), det = TRUE)
        }
        if (!is.null(V$crash.parms)) {
            cat("varcov.spatial: improper matrix for following the given parameters:")
            print(V$crash.parms)
            stop()
        }
        ivx <- solve(V$varcov, xmat[[i]])
        xivx <- crossprod(ivx, xmat[[i]])
        betahat <- .solve.geoR(xivx, crossprod(ivx, data[[i]]))
        res <- data[[i]]
        ssres <- drop(crossprod(res, solve(V$varcov, res)))
        if (method.lik == "ML") {
            negloglik <- (n/2) * (log(2 * pi)) + V$log.det.to.half + 
                0.5 * ssres
        }
        if (method.lik == "RML") {
            choldet <- sum(log(diag(chol(xivx))))
            negloglik <- V$log.det.to.half + 0.5 * ssres + choldet
            xx.eigen <- eigen(crossprod(xmat[[i]]), symmetric = TRUE, 
                only.values = TRUE)
            negloglik <- negloglik + ((n - beta.size)/2) * (log(2 * 
                pi)) - 0.5 * sum(log(xx.eigen$values))
        }
        sumnegloglik <- sumnegloglik + negloglik
    }
    sumnegloglik <- sumnegloglik - log.jacobian
    if (sumnegloglik > (.Machine$double.xmax^0.5)) 
        sumnegloglik <- .Machine$double.xmax^0.5
    return(as.vector(-sumnegloglik))
}





#############################
# no Vecchia parallel MCL_0 #
#############################
MCL0_p_exp <- function(grid, Z, center, border, core){
 
  N <- nrow(grid)
  K <- nrow(center)
  Z <- order.fun(Z, grid, center)                  
  grid <- order.fun(grid, grid, center)             
  grid_region <- region.fun(grid, center)                
  Q <- matrix(,2,K)
  q <- sapply(1:K, function(k){sum(grid_region==k)})
  Q[1,] <- cumsum(q)-q+1
  Q[2,] <- cumsum(q)

  ####### initial ######
  tau <- 0.1
  sigmasqvector <- matrix(c(1,1,0.25,0.25),K,1)
  alphavector <- matrix(c(16,4,16,4),K,1)
  nuvector <- matrix(rep(0.5,K), K, 1)

  ###### Step1 ######
  tt <- proc.time()
  export_vars <- c("loglik.GRF2")
  mclu <- makeCluster(core)
  registerDoParallel(mclu)    
  R1 <- foreach(v=1:K, .combine="cbind", .packages = c("stats", "geoR"), .export = export_vars) %dopar% {         
          like <- function(par){
            loglik.GRF2(coords = matrix(grid[Q[1,v]:Q[2,v],],q[v],ncol(grid)), data = matrix(Z[Q[1,v]:Q[2,v]],q[v],1), 
            obj.model = NULL, cov.model = "exp", cov.pars=c(par[1],1/par[2]), nugget = par[3], 
            kappa = 0.5, lambda = 1, psiR = 1, psiA = 0, trend = "cte", 
            method.lik = "ML", compute.dists = TRUE, realisations = NULL) }
          E <- optim(par = c(sigmasqvector[v], alphavector[v], tau), like, lower = c(rep(1e-06,3)), upper = c(rep(Inf,3)), 
                 method = "L-BFGS-B", control = list(fnscale=-1))
          return(E$par)
          }
  stopCluster(mclu) 
  tau_hat <- sum(R1[3,]*q/N)
  t1 <- proc.time()-tt
  #print(t1)

  ###### Step2 ######
  tt <- proc.time()
  export_vars <- c("loglik.GRF2")
  mclu <- makeCluster(core)
  registerDoParallel(mclu)    
  R2 <- foreach(v=1:K, .combine="cbind", .packages = c("stats", "geoR"), .export = export_vars) %dopar% { 
         like <- function(par){
           loglik.GRF2(coords = matrix(grid[Q[1,v]:Q[2,v],],q[v],ncol(grid)), data = matrix(Z[Q[1,v]:Q[2,v]],q[v],1), 
           obj.model = NULL, cov.model = "exp", cov.pars=c(par[1],1/par[2]), nugget = tau_hat, 
           kappa = 0.5, lambda = 1, psiR = 1, psiA = 0, trend = "cte", 
           method.lik = "ML", compute.dists = TRUE, realisations = NULL) }
         E <- optim(par =c(R1[1,v]+R1[3,v]-tau_hat, R1[2,v]), like, lower = c(rep(1e-06,2)), upper = c(rep(Inf,2)), 
           method = "L-BFGS-B", control = list(fnscale=-1))
         return(E$par)
         }
  stopCluster(mclu) 
  sigmasqvector_hat <- matrix(R2[1,],K,1)
  alphavector_hat <- matrix(R2[2,],K,1)
  nuvector_hat <- nuvector
  t2 <- proc.time()-tt
  #print(t2)
  
  ###### Step3 ######
  tt <- proc.time()
  first1000 <- list()  
  for(h in 1:K){
    if(q[h]>1000){
       first1000[[h]] <- sample(Q[1,h]:Q[2,h], size =1000, replace = F)
    }else{
       first1000[[h]] <- Q[1,h]:Q[2,h]
    }
  }

  cl <- function (par){
          export_vars <- c("alpha.fun", "nu.fun", "rho.fun", "region.fun", "neighbor.fun", "d.fun", "weight.fun1", "weight.fun2",
                           "cov.fun.4", "tau_hat", "sigmasqvector_hat", "alphavector_hat", "nuvector_hat", "center", "border",
                           "first1000", "grid", "Z")
          mclu <- makeCluster(core)
          registerDoParallel(mclu)   
          xx <- foreach(v=1:K, .combine="cbind", .packages = c("SpatialTools", "geoR", "expint", "deldir", "mvnfast"), .export = export_vars) %dopar% { 
                  L <- length(first1000[[v]])
                  error <- diag(tau_hat, L)
                  datacovariance <- cov.fun.4(matrix(grid[first1000[[v]],],L,ncol(grid)), 
                                              center, border, sigmasqvector_hat, alphavector_hat, nuvector_hat, par[1]) + error    
                  dmvn(t(Z[first1000[[v]]]), rep(0,L), datacovariance, log = T)  }
          stopCluster(mclu)
          compo <- sum( xx ) 
    
          #print(-compo)
          return(-compo)
          stopCluster(mclu) 
        }     
  R3 <- optimize(cl, c(0,0.5), tol=1e-06) 
  a_hat <- R3$minimum
  cl_value <- -R3$objective
  t3 <- proc.time()-tt
  #print(t3)

  return( list( center = center,
                border = border,
                Z = Z,
                grid = grid,
                beta0 = 0,
                sigmasqvector_hat = sigmasqvector_hat,
                alphavector_hat = alphavector_hat,
                nuvector_hat = nuvector_hat,
                a_hat = a_hat,
                tau_hat = tau_hat,
                cl = cl_value,
                time = rbind(t1,t2,t3) ) )
}


###################
# no Vecchia MCL_0#
###################
MCL0_exp <- function(grid, Z, center, border){
 
  N <- nrow(grid)
  K <- nrow(center)
  Z <- order.fun(Z, grid, center)                  
  grid <- order.fun(grid, grid, center)             
  grid_region <- region.fun(grid, center)                
  Q <- matrix(,2,K)
  q <- sapply(1:K, function(k){sum(grid_region==k)})
  Q[1,] <- cumsum(q)-q+1
  Q[2,] <- cumsum(q)

  ####### initial ######
  tau <- 0.1
  sigmasqvector <- matrix(c(1,1,0.25,0.25),K,1)
  alphavector <- matrix(c(16,4,16,4),K,1)
  nuvector <- matrix(rep(0.5,K), K, 1)

  ###### Step1 ######
  tt <- proc.time()   
  R1 <- sapply(1:K, function(v){        
          like <- function(par){
            loglik.GRF2(coords = matrix(grid[Q[1,v]:Q[2,v],],q[v],ncol(grid)), data = matrix(Z[Q[1,v]:Q[2,v]],q[v],1), 
            obj.model = NULL, cov.model = "exp", cov.pars=c(par[1],1/par[2]), nugget = par[3], 
            kappa = 0.5, lambda = 1, psiR = 1, psiA = 0, trend = "cte", 
            method.lik = "ML", compute.dists = TRUE, realisations = NULL) }
          E <- optim(par = c(sigmasqvector[v], alphavector[v], tau), like, lower = c(rep(1e-06,3)), upper = c(rep(Inf,3)), 
                 method = "L-BFGS-B", control = list(fnscale=-1))
          return(E$par)
          })
  tau_hat <- sum(R1[3,]*q/N)
  t1 <- proc.time()-tt
  #print(t1)

  ###### Step2 ######
  tt <- proc.time()    
  R2 <- sapply(1:K, function(v){ 
         like <- function(par){
           loglik.GRF2(coords = matrix(grid[Q[1,v]:Q[2,v],],q[v],ncol(grid)), data = matrix(Z[Q[1,v]:Q[2,v]],q[v],1), 
           obj.model = NULL, cov.model = "exp", cov.pars=c(par[1],1/par[2]), nugget = tau_hat, 
           kappa = 0.5, lambda = 1, psiR = 1, psiA = 0, trend = "cte", 
           method.lik = "ML", compute.dists = TRUE, realisations = NULL) }
         E <- optim(par =c(R1[1,v]+R1[3,v]-tau_hat, R1[2,v]), like, lower = c(rep(1e-06,2)), upper = c(rep(Inf,2)), 
           method = "L-BFGS-B", control = list(fnscale=-1))
         return(E$par)
         })
  sigmasqvector_hat <- matrix(R2[1,],K,1)
  alphavector_hat <- matrix(R2[2,],K,1)
  nuvector_hat <- nuvector
  t2 <- proc.time()-tt
  #print(t2)
  
  ###### Step3 ######
  tt <- proc.time()
  first1000 <- list()  
  for(h in 1:K){
    if(q[h]>1000){
       first1000[[h]] <- sample(Q[1,h]:Q[2,h], size =1000, replace = F)
    }else{
       first1000[[h]] <- Q[1,h]:Q[2,h]
    }
  }

  cl <- function (par){ 
          xx <- sapply(1:K, function(v){ 
                  L <- length(first1000[[v]])
                  error <- diag(tau_hat, L)
                  datacovariance <- cov.fun.4(matrix(grid[first1000[[v]],],L,ncol(grid)), 
                                              center, border, sigmasqvector_hat, alphavector_hat, nuvector_hat, par[1]) + error    
                  dmvn(t(Z[first1000[[v]]]), rep(0,L), datacovariance, log = T)  })
          compo <- sum( xx ) 
          return(-compo)
        }     
  R3 <- optimize(cl, c(0,0.5), tol=1e-06) 
  a_hat <- R3$minimum
  cl_value <- -R3$objective
  t3 <- proc.time()-tt
  #print(t3)

  return( list( center = center,
                border = border,
                Z = Z,
                grid = grid,
                beta0 = 0,
                sigmasqvector_hat = sigmasqvector_hat,
                alphavector_hat = alphavector_hat,
                nuvector_hat = nuvector_hat,
                a_hat = a_hat,
                tau_hat = tau_hat,
                cl = cl_value,
                time = rbind(t1,t2,t3) ) )
}


###########################
# no Vecchia parallel MCL #
###########################
MCL_p_exp <- function(grid, Z, center, border, core){
 
  N <- nrow(grid)
  K <- nrow(center)
  Z <- order.fun(Z, grid, center)                  
  grid <- order.fun(grid, grid, center)             
  grid_region <- region.fun(grid, center)                
  Q <- matrix(,2,K)
  q <- sapply(1:K, function(k){sum(grid_region==k)})
  Q[1,] <- cumsum(q)-q+1
  Q[2,] <- cumsum(q)

  ####### initial ######
  tau <- 0.1
  sigmasqvector <- matrix(c(1,1,0.25,0.25),K,1)
  alphavector <- matrix(c(16,4,16,4),K,1)
  nuvector <- matrix(rep(0.5,K), K, 1)

  ###### Step1 ######
  tt <- proc.time()
  export_vars <- c("loglik.GRF2")
  mclu <- makeCluster(core)
  registerDoParallel(mclu)    
  R1 <- foreach(v=1:K, .combine="cbind", .packages = c("stats", "geoR"), .export = export_vars) %dopar% {         
          like <- function(par){
            loglik.GRF2(coords = matrix(grid[Q[1,v]:Q[2,v],],q[v],ncol(grid)), data = matrix(Z[Q[1,v]:Q[2,v]],q[v],1), 
            obj.model = NULL, cov.model = "exp", cov.pars=c(par[1],1/par[2]), nugget = par[3], 
            kappa = 0.5, lambda = 1, psiR = 1, psiA = 0, trend = "cte", 
            method.lik = "ML", compute.dists = TRUE, realisations = NULL) }
          E <- optim(par = c(sigmasqvector[v], alphavector[v], tau), like, lower = c(rep(1e-06,3)), upper = c(rep(Inf,3)), 
                 method = "L-BFGS-B", control = list(fnscale=-1))
          return(E$par)
          }
  stopCluster(mclu) 
  tau_hat <- sum(R1[3,]*q/N)
  t1 <- proc.time()-tt
  #print(t1)

  ###### Step2 ######
  tt <- proc.time()
  export_vars <- c("loglik.GRF2")
  mclu <- makeCluster(core)
  registerDoParallel(mclu)    
  R2 <- foreach(v=1:K, .combine="cbind", .packages = c("stats", "geoR"), .export = export_vars) %dopar% { 
         like <- function(par){
           loglik.GRF2(coords = matrix(grid[Q[1,v]:Q[2,v],],q[v],ncol(grid)), data = matrix(Z[Q[1,v]:Q[2,v]],q[v],1), 
           obj.model = NULL, cov.model = "exp", cov.pars=c(par[1],1/par[2]), nugget = tau_hat, 
           kappa = 0.5, lambda = 1, psiR = 1, psiA = 0, trend = "cte", 
           method.lik = "ML", compute.dists = TRUE, realisations = NULL) }
         E <- optim(par =c(R1[1,v]+R1[3,v]-tau_hat, R1[2,v]), like, lower = c(rep(1e-06,2)), upper = c(rep(Inf,2)), 
           method = "L-BFGS-B", control = list(fnscale=-1))
         return(E$par)
         }
  stopCluster(mclu) 
  sigmasqvector_hat <- matrix(R2[1,],K,1)
  alphavector_hat <- matrix(R2[2,],K,1)
  nuvector_hat <- nuvector
  t2 <- proc.time()-tt
  #print(t2)
  
  ###### Step3 ######
  tt <- proc.time()
  first1000 <- list()  
  for(h in 1:K){
    if(q[h]>1000){
       first1000[[h]] <- sample(Q[1,h]:Q[2,h], size =1000, replace = F)
    }else{
       first1000[[h]] <- Q[1,h]:Q[2,h]
    }
  }

  cl <- function (par){
          export_vars <- c("alpha.fun", "nu.fun", "rho.fun", "region.fun", "neighbor.fun", "d.fun", "weight.fun1", "weight.fun2",
                           "cov.fun.4", "tau_hat", "sigmasqvector_hat", "alphavector_hat", "nuvector_hat", "center", "border",
                           "first1000", "grid", "Z")
          mclu <- makeCluster(core)
          registerDoParallel(mclu)   
          xx <- foreach(v=1:K, .combine="cbind", .packages = c("SpatialTools", "geoR", "expint", "deldir", "mvnfast"), .export = export_vars) %dopar% { 
                  L <- length(first1000[[v]])
                  error <- diag(tau_hat, L)
                  datacovariance <- cov.fun.4(matrix(grid[first1000[[v]],],L,ncol(grid)), 
                                              center, border, sigmasqvector_hat, alphavector_hat, nuvector_hat, par[1]) + error    
                  dmvn(t(Z[first1000[[v]]]), rep(0,L), datacovariance, log = T)  }
          stopCluster(mclu)
          compo <- sum( xx ) 
    
          #print(-compo)
          return(-compo)
          stopCluster(mclu) 
        }     
  R3 <- optimize(cl, c(0,0.5), tol=1e-06) 
  a_hat <- R3$minimum
  t3 <- proc.time()-tt
  #print(t3)

  ###### Step4 ######
  tt <- proc.time()
  lambdavector_hat <- matrix(sigmasqvector_hat*alphavector_hat^(2*nuvector_hat),K,1)
  cl <- function (par){
          export_vars <- c("alpha.fun", "nu.fun", "rho.fun", "region.fun", "neighbor.fun", "d.fun", "weight.fun1", "weight.fun2",
                           "cov.fun.4", "tau_hat", "sigmasqvector_hat", "alphavector_hat", "nuvector_hat", "a_hat", "center", "border",
                           "first1000", "grid", "Z", "K", "lambdavector_hat")
          mclu <- makeCluster(core)
          registerDoParallel(mclu)   
          alpha <- matrix((lambdavector_hat/par[3:(K+2)])^(1/(2*nuvector_hat)),K,1)
          xx <- foreach(v=1:K, .combine="cbind", .packages = c("SpatialTools", "geoR", "expint", "deldir", "mvnfast"), .export = export_vars) %dopar% { 
                  L <- length(first1000[[v]])
                  error <- diag(par[1], L)
                  datacovariance <- cov.fun.4(matrix(grid[first1000[[v]],],L,ncol(grid)), 
                                              center, border, par[3:(K+2)], alpha, nuvector_hat, par[2]) + error    
                  dmvn(t(Z[first1000[[v]]]), rep(0,L), datacovariance, log = T)  }
          stopCluster(mclu)
          compo <- sum( xx ) 
    
          #print(compo)
          return(compo)
          stopCluster(mclu) 
        }     
  R4 <- optim(par = c(tau_hat, a_hat, sigmasqvector_hat), cl, lower = c(1e-06, 0, rep(1e-06,K)), upper = c(Inf, 0.5, rep(Inf,K)), 
            method = "L-BFGS-B", control = list(fnscale=-1,factr=1e12)) 

  sigmasqvector_hat <- matrix(R4$par[3:(K+2)],K,1)
  alphavector_hat <- matrix((lambdavector_hat/R4$par[3:(K+2)])^(1/(2*nuvector_hat)),K,1)
  a_hat <- R4$par[2]
  tau_hat <- R4$par[1]
  cl_value <- R4$value
  t4 <- proc.time()-tt
  #print(t4)

  return( list( center = center,
                border = border,
                Z = Z,
                grid = grid,
                beta0 = 0,
                sigmasqvector_hat = sigmasqvector_hat,
                alphavector_hat = alphavector_hat,
                nuvector_hat = nuvector_hat,
                a_hat = a_hat,
                tau_hat = tau_hat,
                cl = cl_value,
                time = rbind(t1,t2,t3,t4) ) )
}

##################
# no Vecchia MCL #
##################
MCL_exp <- function(grid, Z, center, border){
 
  N <- nrow(grid)
  K <- nrow(center)
  Z <- order.fun(Z, grid, center)                  
  grid <- order.fun(grid, grid, center)             
  grid_region <- region.fun(grid, center)                
  Q <- matrix(,2,K)
  q <- sapply(1:K, function(k){sum(grid_region==k)})
  Q[1,] <- cumsum(q)-q+1
  Q[2,] <- cumsum(q)

  ####### initial ######
  tau <- 0.1
  sigmasqvector <- matrix(c(1,1,0.25,0.25),K,1)
  alphavector <- matrix(c(16,4,16,4),K,1)
  nuvector <- matrix(rep(0.5,K), K, 1)

  ###### Step1 ######
  tt <- proc.time()   
  R1 <- sapply(1:K, function(v){         
          like <- function(par){
            loglik.GRF2(coords = matrix(grid[Q[1,v]:Q[2,v],],q[v],ncol(grid)), data = matrix(Z[Q[1,v]:Q[2,v]],q[v],1), 
            obj.model = NULL, cov.model = "exp", cov.pars=c(par[1],1/par[2]), nugget = par[3], 
            kappa = 0.5, lambda = 1, psiR = 1, psiA = 0, trend = "cte", 
            method.lik = "ML", compute.dists = TRUE, realisations = NULL) }
          E <- optim(par = c(sigmasqvector[v], alphavector[v], tau), like, lower = c(rep(1e-06,3)), upper = c(rep(Inf,3)), 
                 method = "L-BFGS-B", control = list(fnscale=-1))
          return(E$par)
          }) 
  tau_hat <- sum(R1[3,]*q/N)
  t1 <- proc.time()-tt
  #print(t1)

  ###### Step2 ######
  tt <- proc.time()   
  R2 <- sapply(1:K, function(v){ 
         like <- function(par){
           loglik.GRF2(coords = matrix(grid[Q[1,v]:Q[2,v],],q[v],ncol(grid)), data = matrix(Z[Q[1,v]:Q[2,v]],q[v],1), 
           obj.model = NULL, cov.model = "exp", cov.pars=c(par[1],1/par[2]), nugget = tau_hat, 
           kappa = 0.5, lambda = 1, psiR = 1, psiA = 0, trend = "cte", 
           method.lik = "ML", compute.dists = TRUE, realisations = NULL) }
         E <- optim(par =c(R1[1,v]+R1[3,v]-tau_hat, R1[2,v]), like, lower = c(rep(1e-06,2)), upper = c(rep(Inf,2)), 
           method = "L-BFGS-B", control = list(fnscale=-1))
         return(E$par)
         })
  sigmasqvector_hat <- matrix(R2[1,],K,1)
  alphavector_hat <- matrix(R2[2,],K,1)
  nuvector_hat <- nuvector
  t2 <- proc.time()-tt
  #print(t2)
  
  ###### Step3 ######
  tt <- proc.time()
  first1000 <- list()  
  for(h in 1:K){
    if(q[h]>1000){
       first1000[[h]] <- sample(Q[1,h]:Q[2,h], size =1000, replace = F)
    }else{
       first1000[[h]] <- Q[1,h]:Q[2,h]
    }
  }

  cl <- function (par){ 
          xx <- sapply(1:K, function(v){ 
                  L <- length(first1000[[v]])
                  error <- diag(tau_hat, L)
                  datacovariance <- cov.fun.4(matrix(grid[first1000[[v]],],L,ncol(grid)), 
                                              center, border, sigmasqvector_hat, alphavector_hat, nuvector_hat, par[1]) + error    
                  dmvn(t(Z[first1000[[v]]]), rep(0,L), datacovariance, log = T)  })
          compo <- sum( xx ) 
          return(-compo) 
        }     
  R3 <- optimize(cl, c(0,0.5), tol=1e-06) 
  t3 <- proc.time()-tt
  a_hat <- R3$minimum

  ###### Step4 ######
  tt <- proc.time()
  lambdavector_hat <- matrix(sigmasqvector_hat*alphavector_hat^(2*nuvector_hat),K,1)
  cl <- function (par){   
          alpha <- matrix((lambdavector_hat/par[3:(K+2)])^(1/(2*nuvector_hat)),K,1)
          xx <- sapply(1:K, function(v){ 
                  L <- length(first1000[[v]])
                  error <- diag(par[1], L)
                  datacovariance <- cov.fun.4(matrix(grid[first1000[[v]],],L,ncol(grid)), 
                                              center, border, par[3:(K+2)], alpha, nuvector_hat, par[2]) + error    
                  dmvn(t(Z[first1000[[v]]]), rep(0,L), datacovariance, log = T)  })
          compo <- sum( xx ) 
          return(compo)
        }     
  R4 <- optim(par = c(tau_hat, a_hat, sigmasqvector_hat), cl, lower = c(1e-06, 0, rep(1e-06,K)), upper = c(Inf, 0.5, rep(Inf,K)), 
            method = "L-BFGS-B", control = list(fnscale=-1,factr=1e12)) 
  t4 <- proc.time()-tt
  sigmasqvector_hat <- matrix(R4$par[3:(K+2)],K,1)
  alphavector_hat <- matrix((lambdavector_hat/R4$par[3:(K+2)])^(1/(2*nuvector_hat)),K,1)
  a_hat <- R4$par[2]
  tau_hat <- R4$par[1]
  cl_value <- R4$value

  return( list( center = center,
                border = border,
                Z = Z,
                grid = grid,
                beta0 = 0,
                sigmasqvector_hat = sigmasqvector_hat,
                alphavector_hat = alphavector_hat,
                nuvector_hat = nuvector_hat,
                a_hat = a_hat,
                tau_hat = tau_hat,
                cl = cl_value,
                time = rbind(t1,t2,t3,t4) ) )
}



########################
#Vecchia parallel MCL_0#
########################
MCL0_p_exp.2 <- function(grid, Z, center, border, core){
 
  N <- nrow(grid)
  K <- nrow(center)
  Z <- order.fun(Z, grid, center)                  
  grid <- order.fun(grid, grid, center)             
  grid_region <- region.fun(grid, center)                
  Q <- matrix(,2,K)
  q <- sapply(1:K, function(k){sum(grid_region==k)})
  Q[1,] <- cumsum(q)-q+1
  Q[2,] <- cumsum(q)

  ####### initial ######
  tau <- 0.01
  sigmasqvector <- matrix(rep(1,K), K, 1) 
  alphavector <- matrix(rep(2,K), K, 1)
  nuvector <- matrix(rep(0.5,K), K, 1)

  ###### Step1 ######
  tt <- proc.time()
  mclu <- makeCluster(core)
  registerDoParallel(mclu)    
  R1 <- foreach(v=1:K, .combine="cbind", .packages = c("stats", "geoR")) %dopar% {
          locs <- matrix(grid[Q[1,v]:Q[2,v],],q[v],ncol(grid))
          D <- matrix(Z[Q[1,v]:Q[2,v]],q[v],1)
          ord <- GpGp::order_maxmin(locs) # calculate an ordering
          locsord <- locs[ord,] # reorder locations
          Dord <- D[ord,]
          NNar <- GpGp::find_ordered_nn(locsord,50)
          rm(locs)
          rm(D)
          rm(ord)         
          like <- function(par){
            GpGp::vecchia_meanzero_loglik( c(par[1],1/par[2],par[3]/par[1]), "exponential_isotropic", y=Dord, locs=locsord, NNarray=NNar )
          }
          E <- optim(par = c(sigmasqvector[v], alphavector[v], tau), like, lower = c(rep(1e-06,3)), upper = c(rep(Inf,3)), 
                 method = "L-BFGS-B", control = list(fnscale=-1))
          return(E$par)
          }
  stopCluster(mclu) 
  tau_hat <- sum(R1[3,]*q/N)
  ttt <- proc.time()-tt
  print(ttt)

  ###### Step2 ######
  tt <- proc.time()
  mclu <- makeCluster(core)
  registerDoParallel(mclu)    
  R2 <- foreach(v=1:K, .combine="cbind", .packages = c("stats", "geoR")) %dopar% { 
          locs <- matrix(grid[Q[1,v]:Q[2,v],],q[v],ncol(grid))
          D <- matrix(Z[Q[1,v]:Q[2,v]],q[v],1)
          ord <- GpGp::order_maxmin(locs) # calculate an ordering
          locsord <- locs[ord,] # reorder locations
          Dord <- D[ord,]
          NNar <- GpGp::find_ordered_nn(locsord,50)
          rm(locs)
          rm(D)
          rm(ord)         
          like <- function(par){
            GpGp::vecchia_meanzero_loglik( c(par[1],1/par[2],tau_hat/par[1]), "exponential_isotropic", y=Dord, locs=locsord, NNarray=NNar )
          }
          E <- optim(par =c(R1[1,v]+R1[3,v]-tau_hat, R1[2,v]), like, lower = c(rep(1e-06,2)), upper = c(rep(Inf,2)), 
                 method = "L-BFGS-B", control = list(fnscale=-1))
           return(E$par)
        }
  stopCluster(mclu) 
  sigmasqvector_hat <- matrix(R2[1,],K,1)
  alphavector_hat <- matrix(R2[2,],K,1)
  nuvector_hat <- nuvector
  ttt <- proc.time()-tt
  print(ttt)
  
  ###### Step3 ######
  tt <- proc.time()
  first1000 <- list()  
  for(h in 1:K){
    if(q[h]>1000){
       first1000[[h]] <- sample(Q[1,h]:Q[2,h], size =1000, replace = F)
    }else{
       first1000[[h]] <- Q[1,h]:Q[2,h]
    }
  }

   cl <- function (par){
          export_vars <- c("alpha.fun", "nu.fun", "rho.fun", "region.fun", "neighbor.fun", "d.fun", "weight.fun1", "weight.fun2", "cov.fun.4",
                           "tau_hat", "sigmasqvector_hat", "alphavector_hat", "nuvector_hat", "center", "border", "first1000", "grid", "Z")
          mclu <- makeCluster(core)
          registerDoParallel(mclu)   
          xx <- foreach(v=1:K, .combine="cbind", .packages = c("SpatialTools", "geoR", "expint", "deldir", "mvnfast"), .export = export_vars) %dopar% { 
                  L <- length(first1000[[v]])
                  error <- diag(tau_hat, L)
                  datacovariance <- cov.fun.4(matrix(grid[first1000[[v]],],L,ncol(grid)),  
                                              center, border, sigmasqvector_hat, alphavector_hat, nuvector_hat, par[1]) + error    
                  dmvn(t(Z[first1000[[v]]]), rep(0,L), datacovariance, log = T)  }
          stopCluster(mclu)
          compo <- sum( xx ) 

          print(-compo)
          return(-compo)
          stopCluster(mclu) 
        }     
  R3 <- optimize(cl, c(0,0.5), tol=1e-03) 
  a_hat <- R3$minimum
  cl_value <- -R3$objective
  ttt <- proc.time()-tt
  print(ttt)

  return( list( center = center,
                border = border,
                Z = Z,
                grid = grid,
                beta0 = 0,
                sigmasqvector_hat = sigmasqvector_hat,
                alphavector_hat = alphavector_hat,
                nuvector_hat = nuvector_hat,
                a_hat = a_hat,
                tau_hat = tau_hat,
                cl = cl_value ) )
}



###############
#Vecchia MCL_0#
###############
MCL0_exp.2 <- function(grid, Z, center, border){
 
  N <- nrow(grid)
  K <- nrow(center)
  Z <- order.fun(Z, grid, center)                  
  grid <- order.fun(grid, grid, center)             
  grid_region <- region.fun(grid, center)                
  Q <- matrix(,2,K)
  q <- sapply(1:K, function(k){sum(grid_region==k)})
  Q[1,] <- cumsum(q)-q+1
  Q[2,] <- cumsum(q)

  ####### initial ######
  tau <- 0.01
  sigmasqvector <- matrix(rep(1,K), K, 1) 
  alphavector <- matrix(rep(2,K), K, 1)
  nuvector <- matrix(rep(0.5,K), K, 1)

  ###### Step1 ######
  set.seed(0)
  tt <- proc.time()    
  R1 <- sapply(1:K, function(v){
          locs <- matrix(grid[Q[1,v]:Q[2,v],],q[v],ncol(grid))
          D <- matrix(Z[Q[1,v]:Q[2,v]],q[v],1)
          ord <- GpGp::order_maxmin(locs) # calculate an ordering
          locsord <- locs[ord,] # reorder locations
          Dord <- D[ord,]
          NNar <- GpGp::find_ordered_nn(locsord,50)
          rm(locs)
          rm(D)
          rm(ord)         
          like <- function(par){
            GpGp::vecchia_meanzero_loglik( c(par[1],1/par[2],par[3]/par[1]), "exponential_isotropic", y=Dord, locs=locsord, NNarray=NNar )
            }
          E <- optim(par = c(sigmasqvector[v], alphavector[v], tau), like, lower = c(rep(1e-06,3)), upper = c(rep(Inf,3)), 
                 method = "L-BFGS-B", control = list(fnscale=-1))
          return(E$par)
        }) 
  ttt <- proc.time()-tt
  print(ttt)
  tau_hat <- sum(R1[3,]*q/N)

  ###### Step2 ######
  tt <- proc.time()   
  R2 <- sapply(1:K, function(v){ 
          locs <- matrix(grid[Q[1,v]:Q[2,v],],q[v],ncol(grid))
          D <- matrix(Z[Q[1,v]:Q[2,v]],q[v],1)
          ord <- GpGp::order_maxmin(locs) # calculate an ordering
          locsord <- locs[ord,] # reorder locations
          Dord <- D[ord,]
          NNar <- GpGp::find_ordered_nn(locsord,50)
          rm(locs)
          rm(D)
          rm(ord)         
          like <- function(par){
            GpGp::vecchia_meanzero_loglik( c(par[1],1/par[2],tau_hat/par[1]), "exponential_isotropic", y=Dord, locs=locsord, NNarray=NNar )
            }
          E <- optim(par =c(R1[1,v]+R1[3,v]-tau_hat, R1[2,v]), like, lower = c(rep(1e-06,2)), upper = c(rep(Inf,2)), 
                 method = "L-BFGS-B", control = list(fnscale=-1))
          return(E$par)
        })
  ttt <- proc.time()-tt
  print(ttt) 
  sigmasqvector_hat <- matrix(R2[1,],K,1)
  alphavector_hat <- matrix(R2[2,],K,1)
  nuvector_hat <- nuvector
  
  ###### Step3 ######
  tt <- proc.time()
  first1000 <- list()  
  for(h in 1:K){
    if(q[h]>1000){
       first1000[[h]] <- sample(Q[1,h]:Q[2,h], size =1000, replace = F)
    }else{
       first1000[[h]] <- Q[1,h]:Q[2,h]
    }
  }

   cl <- function (par){  
          xx <- sapply(1:K, function(v){ 
                  L <- length(first1000[[v]])
                  error <- diag(tau_hat, L)
                  datacovariance <- cov.fun.4(matrix(grid[first1000[[v]],],L,ncol(grid)),  
                                              center, border, sigmasqvector_hat, alphavector_hat, nuvector_hat, par[1]) + error    
                  dmvn(t(Z[first1000[[v]]]), rep(0,L), datacovariance, log = T)  })
          compo <- sum( xx )
          print(-compo) 
          return(-compo)
        }     
  R3 <- optimize(cl, c(0,0.5), tol=1e-03) 
  ttt <- proc.time()-tt
  print(ttt)
  a_hat <- R3$minimum
  cl_value <- -R3$objective

  return( list( center = center,
                border = border,
                Z = Z,
                grid = grid,
                beta0 = 0,
                sigmasqvector_hat = sigmasqvector_hat,
                alphavector_hat = alphavector_hat,
                nuvector_hat = nuvector_hat,
                a_hat = a_hat,
                tau_hat = tau_hat,
                cl = cl_value ) )
}



 




#parallel DK

DK_p <- function(E, grid = E$grid, Z = E$Z, center_p, grids, core){
  
  Z_p <- order.fun(Z, grid, center_p)                  
  grid_p <- order.fun(grid, grid, center_p) 
  grid_p_region <- region.fun(grid_p, center_p)
  M <- nrow(center_p)                       
  Qmatrix_p <- matrix(,2,M)
  qq <- sapply(1:M, function(k){sum(grid_p_region==k)})
  Qmatrix_p[1,] <- cumsum(qq)-qq+1
  Qmatrix_p[2,] <- cumsum(qq)

  tt <- proc.time()
  export_vars <- c("alpha.fun", "nu.fun", "rho.fun", "region.fun", "neighbor.fun", "d.fun", "weight.fun1", "weight.fun2", "cov.fun.4")
  mclu <- makeCluster(core)
  registerDoParallel(mclu)   
  Sigma_yii <- foreach(k=1:M, .packages = c("SpatialTools","geoR", "expint", "deldir"), .export = export_vars)%dopar%{
                 error <- diag(E$tau_hat, qq[k])
                 covariance <- cov.fun.4(matrix(grid_p[Qmatrix_p[1,k]:Qmatrix_p[2,k],], qq[k], ncol(grid_p)),
                                         E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat)+ error                        
               } 
  stopCluster(mclu)
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()
  Sigma_yij <- list()
  for(j in 1:(M-1)){
    export_vars <- c("alpha.fun", "nu.fun", "rho.fun", "region.fun", "neighbor.fun", "d.fun", "weight.fun1", "weight.fun2", "cov.fun.3")
    mclu <- makeCluster(min(core,M-j))
    registerDoParallel(mclu)   
    Sigma_yij[[j]] <- foreach(i=(j+1):M, .packages = c("SpatialTools","geoR", "expint", "deldir"), .export = export_vars)%dopar%{
                          covariance <- cov.fun.3(matrix(grid_p[Qmatrix_p[1,j]:Qmatrix_p[2,j],],(qq[j]),ncol(grid_p)),
                                                  matrix(grid_p[Qmatrix_p[1,i]:Qmatrix_p[2,i],],(qq[i]),ncol(grid_p)),
                                                  E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat) } 
    stopCluster(mclu)
  }
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()
  mclu <- makeCluster(core)
  registerDoParallel(mclu)   
  Sigma_yii_inverse <- foreach(i=1:M, .packages = c("SpatialTools","geoR"))%dopar%{
                         chol2inv( chol(Sigma_yii[[i]]) ) 
                         }   
  stopCluster(mclu)
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()
  predictsize <- nrow(grids)
  grids_p_region <- region.fun(grids, center_p)
  grids_label <- sapply(1:predictsize, function(i){ sum( grids_p_region[1:i]==grids_p_region[i] ) })
  pcov <- list()
  for(u in 1:M){
      #print(u)
      o <- which(grids_p_region==u)
      if(length(o)==0){
        pcov[[u]] <- 0
      }else{
        export_vars <- c("alpha.fun", "nu.fun", "rho.fun", "region.fun", "neighbor.fun", "d.fun", "weight.fun1", "weight.fun2", "cov.fun.3")
        mclu <- makeCluster(core)
        registerDoParallel(mclu)   
        pcov[[u]] <- foreach(i=1:M, .packages = c("SpatialTools","geoR", "expint", "deldir"), .export = export_vars)%dopar%{
                       cov.fun.3( matrix(grids[o,],length(o),ncol(grids)), 
                                  matrix(grid_p[Qmatrix_p[1,i]:Qmatrix_p[2,i],],qq[i],ncol(grids)), 
                                  E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat)
                     }   
        stopCluster(mclu)                    
      }   
  }           
  ttt <- proc.time()-tt 
  #print(ttt)
  
  tt <- proc.time()
  export_vars <- c("alpha.fun", "nu.fun", "rho.fun", "region.fun", "neighbor.fun", "d.fun", "weight.fun1", "weight.fun2", "cov.fun.4")
  mclu <- makeCluster(core)
  registerDoParallel(mclu)   
  R_pred <- foreach( k=1:predictsize, .combine="cbind", .packages = c("SpatialTools","geoR", "expint", "deldir"), .export = export_vars) %dopar% {
    u <- grids_p_region[k]
    bi <- lapply(1:M, function (i) { matrix(pcov[[ u ]][[ i ]][grids_label[k],] %*% Sigma_yii_inverse[[ i ]], qq[i], 1) } )
    hat_z_local <- matrix( sapply(1:M, function (i) { crossprod(bi[[i]] , matrix(Z_p[Qmatrix_p[1,i]:Qmatrix_p[2,i]], qq[i], 1)) }), M, 1)
    di <- matrix( sapply(1:M, function (i) { tcrossprod( t(bi[[i]]) , pcov[[ u ]][[ i ]][grids_label[k],] ) }), 1, M )
    MM <- matrix(, M, M)
    for (i in 1:M) {
      for (j in 1:M) {
        if(i==j){ 
           MM[i,j] <- crossprod( bi[[i]] , Sigma_yii[[i]] %*% bi[[j]] )
        }else if (i>j){
           MM[i,j] <- crossprod( bi[[i]] , t(Sigma_yij[[j]][[i-j]]) %*% bi[[j]] )
        }else if (i<j){
           MM[i,j] <- crossprod( bi[[i]] , Sigma_yij[[i]][[j-i]] %*% bi[[j]] )
        }  
      }
    }

    pred <- di %*% chol2inv(chol(MM)) %*% hat_z_local
    va <- cov.fun.4(matrix(grids[k,],1,ncol(grids)),
                    E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat) 
    pred_sd <- sqrt(va - di %*% chol2inv(chol(MM)) %*% t(di))
    return( c(pred, pred_sd) )
  }                       
  stopCluster(mclu)
  ttt <- proc.time()-tt
  #print(ttt)
  return( list(grids = grids, predictor = matrix(R_pred[1,],predictsize,1), sd = matrix(R_pred[2,],predictsize,1) ))
}



DK <- function(E, grid = E$grid, Z = E$Z, center_p, grids){
  
  Z_p <- order.fun(Z, grid, center_p)                  
  grid_p <- order.fun(grid, grid, center_p) 
  grid_p_region <- region.fun(grid_p, center_p)
  M <- nrow(center_p)                       
  Qmatrix_p <- matrix(,2,M)
  qq <- sapply(1:M, function(k){sum(grid_p_region==k)})
  Qmatrix_p[1,] <- cumsum(qq)-qq+1
  Qmatrix_p[2,] <- cumsum(qq)

  tt <- proc.time()   
  Sigma_yii <- lapply(1:M, function(k){
                 error <- diag(E$tau_hat, qq[k])
                 covariance <- cov.fun.4(matrix(grid_p[Qmatrix_p[1,k]:Qmatrix_p[2,k],], qq[k], ncol(grid_p)),
                                         E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat)+ error                        
               }) 
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()
  Sigma_yij <- list()
  for(j in 1:(M-1)){   
    Sigma_yij[[j]] <- lapply((j+1):M, function(i){
                          covariance <- cov.fun.3(matrix(grid_p[Qmatrix_p[1,j]:Qmatrix_p[2,j],],(qq[j]),ncol(grid_p)),
                                                  matrix(grid_p[Qmatrix_p[1,i]:Qmatrix_p[2,i],],(qq[i]),ncol(grid_p)),
                                                  E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat) }) 
  }
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()  
  Sigma_yii_inverse <- lapply(1:M, function(i){
                         chol2inv( chol(Sigma_yii[[i]]) ) 
                       })  
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()
  predictsize <- nrow(grids)
  grids_p_region <- region.fun(grids, center_p)
  grids_label <- sapply(1:predictsize, function(i){ sum( grids_p_region[1:i]==grids_p_region[i] ) })
  pcov <- list()
  for(u in 1:M){
      #print(u)
      o <- which(grids_p_region==u)
      if(length(o)==0){
        pcov[[u]] <- 0
      }else{
        pcov[[u]] <- lapply(1:M, function(i){
                       cov.fun.3( matrix(grids[o,],length(o),ncol(grids)), 
                                  matrix(grid_p[Qmatrix_p[1,i]:Qmatrix_p[2,i],],qq[i],ncol(grids)), 
                                  E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat)
                     })                      
      }   
  }           
  ttt <- proc.time()-tt 
  #print(ttt)
  
  tt <- proc.time()  
  R_pred <- sapply(1:predictsize, function(k){
    u <- grids_p_region[k]
    bi <- lapply(1:M, function (i) { matrix(pcov[[ u ]][[ i ]][grids_label[k],] %*% Sigma_yii_inverse[[ i ]], qq[i], 1) } )
    hat_z_local <- matrix( sapply(1:M, function (i) { crossprod(bi[[i]] , matrix(Z_p[Qmatrix_p[1,i]:Qmatrix_p[2,i]], qq[i], 1)) }), M, 1)
    di <- matrix( sapply(1:M, function (i) { tcrossprod( t(bi[[i]]) , pcov[[ u ]][[ i ]][grids_label[k],] ) }), 1, M )
    MM <- matrix(, M, M)
    for (i in 1:M) {
      for (j in 1:M) {
        if(i==j){ 
           MM[i,j] <- crossprod( bi[[i]] , Sigma_yii[[i]] %*% bi[[j]] )
        }else if (i>j){
           MM[i,j] <- crossprod( bi[[i]] , t(Sigma_yij[[j]][[i-j]]) %*% bi[[j]] )
        }else if (i<j){
           MM[i,j] <- crossprod( bi[[i]] , Sigma_yij[[i]][[j-i]] %*% bi[[j]] )
        }  
      }
    }

    pred <- di %*% chol2inv(chol(MM)) %*% hat_z_local
    va <- cov.fun.4(matrix(grids[k,],1,ncol(grids)),
                    E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat) 
    pred_sd <- sqrt(va - di %*% chol2inv(chol(MM)) %*% t(di))
    return( c(pred, pred_sd) )
  })                       
  ttt <- proc.time()-tt
  #print(ttt)
  return( list(grids = grids, predictor = matrix(R_pred[1,],predictsize,1), sd = matrix(R_pred[2,],predictsize,1) ))
}


#not calculate pred_sd

DK.sim <- function(E, grid = E$grid, Z = E$Z, center_p, grids){
  
  Z_p <- order.fun(Z, grid, center_p)                  
  grid_p <- order.fun(grid, grid, center_p) 
  grid_p_region <- region.fun(grid_p, center_p)
  M <- nrow(center_p)                       
  Qmatrix_p <- matrix(,2,M)
  qq <- sapply(1:M, function(k){sum(grid_p_region==k)})
  Qmatrix_p[1,] <- cumsum(qq)-qq+1
  Qmatrix_p[2,] <- cumsum(qq)

  tt <- proc.time()   
  Sigma_yii <- lapply(1:M, function(k){
                 error <- diag(E$tau_hat, qq[k])
                 covariance <- cov.fun.4(matrix(grid_p[Qmatrix_p[1,k]:Qmatrix_p[2,k],], qq[k], ncol(grid_p)),
                                         E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat)+ error                        
               }) 
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()
  Sigma_yij <- list()
  for(j in 1:(M-1)){   
    Sigma_yij[[j]] <- lapply((j+1):M, function(i){
                          covariance <- cov.fun.3(matrix(grid_p[Qmatrix_p[1,j]:Qmatrix_p[2,j],],(qq[j]),ncol(grid_p)),
                                                  matrix(grid_p[Qmatrix_p[1,i]:Qmatrix_p[2,i],],(qq[i]),ncol(grid_p)),
                                                  E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat) }) 
  }
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()  
  Sigma_yii_inverse <- lapply(1:M, function(i){
                         chol2inv( chol(Sigma_yii[[i]]) ) 
                       })  
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()
  predictsize <- nrow(grids)
  grids_p_region <- region.fun(grids, center_p)
  grids_label <- sapply(1:predictsize, function(i){ sum( grids_p_region[1:i]==grids_p_region[i] ) })
  pcov <- list()
  for(u in 1:M){
      #print(u)
      o <- which(grids_p_region==u)
      if(length(o)==0){
        pcov[[u]] <- 0
      }else{
        pcov[[u]] <- lapply(1:M, function(i){
                       cov.fun.3( matrix(grids[o,],length(o),ncol(grids)), 
                                  matrix(grid_p[Qmatrix_p[1,i]:Qmatrix_p[2,i],],qq[i],ncol(grids)), 
                                  E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat)
                     })                      
      }   
  }           
  ttt <- proc.time()-tt 
  #print(ttt)
  
  tt <- proc.time()  
  R_pred <- sapply(1:predictsize, function(k){
    u <- grids_p_region[k]
    bi <- lapply(1:M, function (i) { matrix(pcov[[ u ]][[ i ]][grids_label[k],] %*% Sigma_yii_inverse[[ i ]], qq[i], 1) } )
    hat_z_local <- matrix( sapply(1:M, function (i) { crossprod(bi[[i]] , matrix(Z_p[Qmatrix_p[1,i]:Qmatrix_p[2,i]], qq[i], 1)) }), M, 1)
    di <- matrix( sapply(1:M, function (i) { tcrossprod( t(bi[[i]]) , pcov[[ u ]][[ i ]][grids_label[k],] ) }), 1, M )
    MM <- matrix(, M, M)
    for (i in 1:M) {
      for (j in 1:M) {
        if(i==j){ 
           MM[i,j] <- crossprod( bi[[i]] , Sigma_yii[[i]] %*% bi[[j]] )
        }else if (i>j){
           MM[i,j] <- crossprod( bi[[i]] , t(Sigma_yij[[j]][[i-j]]) %*% bi[[j]] )
        }else if (i<j){
           MM[i,j] <- crossprod( bi[[i]] , Sigma_yij[[i]][[j-i]] %*% bi[[j]] )
        }  
      }
    }

    pred <- di %*% chol2inv(chol(MM)) %*% hat_z_local
    #va <- cov.fun.4(matrix(grids[k,],1,ncol(grids)),
    #                E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat) 
    #pred_sd <- sqrt(va - di %*% chol2inv(chol(MM)) %*% t(di))
    #return( c(pred, pred_sd) )
  })                       
  ttt <- proc.time()-tt
  #print(ttt)
  return( list(grids = grids, predictor = matrix(R_pred,predictsize,1) ))
}


#parallel LDK

LDK_p <- function(E, grid = E$grid, Z = E$Z, center_p, grids, core){

  Z_p <- order.fun(Z, grid, center_p)                  
  grid_p <- order.fun(grid, grid, center_p) 
  grid_p_region <- region.fun(grid_p, center_p)
  M <- nrow(center_p)                       
  Qmatrix_p <- matrix(,2,M)
  qq <- sapply(1:M, function(k){sum(grid_p_region==k)})
  Qmatrix_p[1,] <- cumsum(qq)-qq+1
  Qmatrix_p[2,] <- cumsum(qq)

  vor_dxy_p <- deldir(center_p[,1], center_p[,2], rw = E$border)
  neighbor_p <- neighbor.fun(center_p, E$border)
  Neighbor_p <- lapply(1:M, function(i){
    A <- which( neighbor_p[,i]!=0 )
    B <- rbind(as.matrix(vor_dxy_p$dirsgs[A,1:2]),as.matrix(vor_dxy_p$dirsgs[A,3:4]))
    D <- dist1(B)
    del <- c()
    for(j in 1:dim(B)[1]){
      same <- which(D[j,]==0)
      del <-c( del,same[which((same-j)>0)] )
    }
    if(length(del)!=0){
        B <- B[-del, ]
    }
    E1 <- dist2(B,as.matrix(vor_dxy_p$dirsgs[,1:2]))
    E2 <- dist2(B,as.matrix(vor_dxy_p$dirsgs[,3:4]))
    F1 <- sapply(1:dim(neighbor_p)[1], function(j){prod(E1[,j])==0})
    F2 <- sapply(1:dim(neighbor_p)[1], function(j){prod(E2[,j])==0})
    relation <- which(F1+F2 != 0)
    G <- sapply(1:M, function(k){sum(neighbor_p[relation,k])})
    return( setdiff(which( G!=0 ),i) )
  })

  tt <- proc.time()
  export_vars <- c("alpha.fun", "nu.fun", "rho.fun", "region.fun", "neighbor.fun", "d.fun", "weight.fun1", "weight.fun2", "cov.fun.4")
  mclu <- makeCluster(core)
  registerDoParallel(mclu)   
  Sigma_yii <- foreach(k=1:M, .packages = c("SpatialTools","geoR", "expint", "deldir"), .export = export_vars)%dopar%{
                 error <- diag(E$tau_hat, qq[k])
                 covariance <- cov.fun.4(matrix(grid_p[Qmatrix_p[1,k]:Qmatrix_p[2,k],], qq[k], ncol(grid_p)),
                                         E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat)+ error                        
               } 
  stopCluster(mclu)
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()
  Sigma_yij <- list()
  for(j in 1:(M-1)){
    N <- length(Neighbor_p[[j]])
    A1=A <- Neighbor_p[[j]]
    for(i in 1:N){
      A <- unique(c(A, Neighbor_p[[A1[i]]])) 
    }
    A <- sort(A) 
    export_vars <- c("alpha.fun", "nu.fun", "rho.fun", "region.fun", "neighbor.fun", "d.fun", "weight.fun1", "weight.fun2", "cov.fun.3")
    mclu <- makeCluster(min(core,M-j))
    registerDoParallel(mclu)   
    Sigma_yij[[j]] <- foreach(i=(j+1):M, .packages = c("SpatialTools","geoR", "expint", "deldir"), .export = export_vars)%dopar%{ 
                        if(sum(i == A)){ 
                          covariance <- cov.fun.3(matrix(grid_p[Qmatrix_p[1,j]:Qmatrix_p[2,j],],(qq[j]),ncol(grid_p)),
                                                  matrix(grid_p[Qmatrix_p[1,i]:Qmatrix_p[2,i],],(qq[i]),ncol(grid_p)),
                                                  E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat)
                        }else{
                          return(0)
                        }                      
                      }
    stopCluster(mclu)
  }
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()
  mclu <- makeCluster(core)
  registerDoParallel(mclu)   
  Sigma_yii_inverse <- foreach(i=1:M, .packages = c("SpatialTools","geoR"))%dopar%{
                         chol2inv( chol(Sigma_yii[[i]]) ) 
                         }   
  stopCluster(mclu)
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()
  predictsize <- nrow(grids)
  grids_p_region <- region.fun(grids, center_p)
  grids_label <- sapply(1:predictsize, function(i){ sum( grids_p_region[1:i]==grids_p_region[i] ) })
  pcov <- list()
  for(u in 1:M){
      #print(u)
      o <- which(grids_p_region==u)
      if(length(o)==0){
        pcov[[u]] <- 0
      }else{
        A <- sort(c(u,Neighbor_p[[u]]))
        L <- length(Neighbor_p[[u]])+1
        export_vars <- c("alpha.fun", "nu.fun", "rho.fun", "region.fun", "neighbor.fun", "d.fun", "weight.fun1", "weight.fun2", "cov.fun.3")
        mclu <- makeCluster(min(core,L))
        registerDoParallel(mclu)   
        pcov[[u]] <- foreach(i=1:L, .packages = c("SpatialTools","geoR", "expint", "deldir"), .export = export_vars)%dopar%{
                       cov.fun.3( matrix(grids[o,],length(o),ncol(grids)), 
                                  matrix(grid_p[Qmatrix_p[1,A[i]]:Qmatrix_p[2,A[i]],],qq[A[i]],ncol(grids)), 
                                  E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat)
                     }   
        stopCluster(mclu)                    
      }   
  }           
  ttt <- proc.time()-tt 
  #print(ttt)

  tt <- proc.time()
  export_vars <- c("alpha.fun", "nu.fun", "rho.fun", "region.fun", "neighbor.fun", "d.fun", "weight.fun1", "weight.fun2", "cov.fun.4")
  mclu <- makeCluster(core)
  registerDoParallel(mclu)   
  R_pred <- foreach( k=1:predictsize, .combine="cbind", .packages = c("SpatialTools","geoR", "expint", "deldir"), .export = export_vars) %dopar% {
    #print(k)
    u <- grids_p_region[k]
    A <- sort(c(u,Neighbor_p[[u]]))
    L <- length(Neighbor_p[[u]])+1                   
    bi <- lapply(1:L, function (i) { matrix(pcov[[ u ]][[ i ]][grids_label[k],] %*% Sigma_yii_inverse[[A[i]]], qq[A[i]], 1) } )
    hat_z_local <- matrix( sapply(1:L, function (i) { crossprod(bi[[i]] , matrix(Z_p[Qmatrix_p[1,A[i]]:Qmatrix_p[2,A[i]]], qq[A[i]], 1)) }), L, 1)
    di <- matrix( sapply(1:L, function (i) { tcrossprod( t(bi[[i]]) , pcov[[ u ]][[ i ]][grids_label[k],] ) }), 1, L )
    MM <- matrix(, L, L)
    for (i in 1:L) {
      for (j in 1:L) {
        if(i==j){ 
           MM[i,j] <- crossprod( bi[[i]] , Sigma_yii[[A[i]]] %*% bi[[j]] )
        }else if (A[i]>A[j]){
           MM[i,j] <- crossprod( bi[[i]] , t(Sigma_yij[[A[j]]][[A[i]-A[j]]]) %*% bi[[j]] )
        }else if (A[i]<A[j]){
           MM[i,j] <- crossprod( bi[[i]] , Sigma_yij[[A[i]]][[A[j]-A[i]]] %*% bi[[j]] )
        }  
      }
    }
    pred <- di %*% chol2inv(chol(MM)) %*% hat_z_local
    va <- cov.fun.4(matrix(grids[k,],1,ncol(grids)),
                    E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat) 
    pred_sd <- sqrt(va - di %*% chol2inv(chol(MM)) %*% t(di))
    return( c(pred, pred_sd) )
  }
  stopCluster(mclu)
  ttt <- proc.time()-tt
  #print(ttt)
  return( list(grids = grids, predictor = matrix(R_pred[1,],predictsize,1), sd = matrix(R_pred[2,],predictsize,1) ))
}


LDK <- function(E, grid = E$grid, Z = E$Z, center_p, grids){

  Z_p <- order.fun(Z, grid, center_p)                  
  grid_p <- order.fun(grid, grid, center_p) 
  grid_p_region <- region.fun(grid_p, center_p)
  M <- nrow(center_p)                       
  Qmatrix_p <- matrix(,2,M)
  qq <- sapply(1:M, function(k){sum(grid_p_region==k)})
  Qmatrix_p[1,] <- cumsum(qq)-qq+1
  Qmatrix_p[2,] <- cumsum(qq)

  vor_dxy_p <- deldir(center_p[,1], center_p[,2], rw = E$border)
  neighbor_p <- neighbor.fun(center_p, E$border)
  Neighbor_p <- lapply(1:M, function(i){
    A <- which( neighbor_p[,i]!=0 )
    B <- rbind(as.matrix(vor_dxy_p$dirsgs[A,1:2]),as.matrix(vor_dxy_p$dirsgs[A,3:4]))
    D <- dist1(B)
    del <- c()
    for(j in 1:dim(B)[1]){
      same <- which(D[j,]==0)
      del <-c( del,same[which((same-j)>0)] )
    }
    if(length(del)!=0){
        B <- B[-del, ]
    }
    E1 <- dist2(B,as.matrix(vor_dxy_p$dirsgs[,1:2]))
    E2 <- dist2(B,as.matrix(vor_dxy_p$dirsgs[,3:4]))
    F1 <- sapply(1:dim(neighbor_p)[1], function(j){prod(E1[,j])==0})
    F2 <- sapply(1:dim(neighbor_p)[1], function(j){prod(E2[,j])==0})
    relation <- which(F1+F2 != 0)
    G <- sapply(1:M, function(k){sum(neighbor_p[relation,k])})
    return( setdiff(which( G!=0 ),i) )
  })

  tt <- proc.time()   
  Sigma_yii <- lapply(1:M, function(k){
                 error <- diag(E$tau_hat, qq[k])
                 covariance <- cov.fun.4(matrix(grid_p[Qmatrix_p[1,k]:Qmatrix_p[2,k],], qq[k], ncol(grid_p)),
                                         E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat)+ error                        
               }) 
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()
  Sigma_yij <- list()
  for(j in 1:(M-1)){
    N <- length(Neighbor_p[[j]])
    A1=A <- Neighbor_p[[j]]
    for(i in 1:N){
      A <- unique(c(A, Neighbor_p[[A1[i]]])) 
    }
    A <- sort(A)   
    Sigma_yij[[j]] <- lapply((j+1):M, function(i){ 
                        if(sum(i == A)){ 
                          covariance <- cov.fun.3(matrix(grid_p[Qmatrix_p[1,j]:Qmatrix_p[2,j],],(qq[j]),ncol(grid_p)),
                                                  matrix(grid_p[Qmatrix_p[1,i]:Qmatrix_p[2,i],],(qq[i]),ncol(grid_p)),
                                                  E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat)
                        }else{
                          return(0)
                        }                      
                      })
  }
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()  
  Sigma_yii_inverse <- lapply(1:M, function(i){
                         chol2inv( chol(Sigma_yii[[i]]) ) 
                       })   
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()
  predictsize <- nrow(grids)
  grids_p_region <- region.fun(grids, center_p)
  grids_label <- sapply(1:predictsize, function(i){ sum( grids_p_region[1:i]==grids_p_region[i] ) })
  pcov <- list()
  for(u in 1:M){
      #print(u)
      o <- which(grids_p_region==u)
      if(length(o)==0){
        pcov[[u]] <- 0
      }else{
        A <- sort(c(u,Neighbor_p[[u]]))
        L <- length(Neighbor_p[[u]])+1   
        pcov[[u]] <- lapply(1:L, function(i){
                       cov.fun.3( matrix(grids[o,],length(o),ncol(grids)), 
                                  matrix(grid_p[Qmatrix_p[1,A[i]]:Qmatrix_p[2,A[i]],],qq[A[i]],ncol(grids)), 
                                  E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat)
                     })                      
      }   
  }           
  ttt <- proc.time()-tt 
  #print(ttt)

  tt <- proc.time()  
  R_pred <- sapply(1:predictsize, function(k){
    u <- grids_p_region[k]
    A <- sort(c(u,Neighbor_p[[u]]))
    L <- length(Neighbor_p[[u]])+1                   
    bi <- lapply(1:L, function (i) { matrix(pcov[[ u ]][[ i ]][grids_label[k],] %*% Sigma_yii_inverse[[A[i]]], qq[A[i]], 1) } )
    hat_z_local <- matrix( sapply(1:L, function (i) { crossprod(bi[[i]] , matrix(Z_p[Qmatrix_p[1,A[i]]:Qmatrix_p[2,A[i]]], qq[A[i]], 1)) }), L, 1)
    di <- matrix( sapply(1:L, function (i) { tcrossprod( t(bi[[i]]) , pcov[[ u ]][[ i ]][grids_label[k],] ) }), 1, L )
    MM <- matrix(, L, L)
    for (i in 1:L) {
      for (j in 1:L) {
        if(i==j){ 
           MM[i,j] <- crossprod( bi[[i]] , Sigma_yii[[A[i]]] %*% bi[[j]] )
        }else if (A[i]>A[j]){
           MM[i,j] <- crossprod( bi[[i]] , t(Sigma_yij[[A[j]]][[A[i]-A[j]]]) %*% bi[[j]] )
        }else if (A[i]<A[j]){
           MM[i,j] <- crossprod( bi[[i]] , Sigma_yij[[A[i]]][[A[j]-A[i]]] %*% bi[[j]] )
        }  
      }
    }
    pred <- di %*% chol2inv(chol(MM)) %*% hat_z_local
    va <- cov.fun.4(matrix(grids[k,],1,ncol(grids)),
                    E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat) 
    pred_sd <- sqrt(va - di %*% chol2inv(chol(MM)) %*% t(di))
    return( c(pred, pred_sd) )
  })
  ttt <- proc.time()-tt
  #print(ttt)
  return( list(grids = grids, predictor = matrix(R_pred[1,],predictsize,1), sd = matrix(R_pred[2,],predictsize,1) ))
}

#not calculate pred_sd

LDK.sim <- function(E, grid = E$grid, Z = E$Z, center_p, grids){

  Z_p <- order.fun(Z, grid, center_p)                  
  grid_p <- order.fun(grid, grid, center_p) 
  grid_p_region <- region.fun(grid_p, center_p)
  M <- nrow(center_p)                       
  Qmatrix_p <- matrix(,2,M)
  qq <- sapply(1:M, function(k){sum(grid_p_region==k)})
  Qmatrix_p[1,] <- cumsum(qq)-qq+1
  Qmatrix_p[2,] <- cumsum(qq)

  vor_dxy_p <- deldir(center_p[,1], center_p[,2], rw = E$border)
  neighbor_p <- neighbor.fun(center_p, E$border)
  Neighbor_p <- lapply(1:M, function(i){
    A <- which( neighbor_p[,i]!=0 )
    B <- rbind(as.matrix(vor_dxy_p$dirsgs[A,1:2]),as.matrix(vor_dxy_p$dirsgs[A,3:4]))
    D <- dist1(B)
    del <- c()
    for(j in 1:dim(B)[1]){
      same <- which(D[j,]==0)
      del <-c( del,same[which((same-j)>0)] )
    }
    if(length(del)!=0){
        B <- B[-del, ]
    }
    E1 <- dist2(B,as.matrix(vor_dxy_p$dirsgs[,1:2]))
    E2 <- dist2(B,as.matrix(vor_dxy_p$dirsgs[,3:4]))
    F1 <- sapply(1:dim(neighbor_p)[1], function(j){prod(E1[,j])==0})
    F2 <- sapply(1:dim(neighbor_p)[1], function(j){prod(E2[,j])==0})
    relation <- which(F1+F2 != 0)
    G <- sapply(1:M, function(k){sum(neighbor_p[relation,k])})
    return( setdiff(which( G!=0 ),i) )
  })

  tt <- proc.time()   
  Sigma_yii <- lapply(1:M, function(k){
                 error <- diag(E$tau_hat, qq[k])
                 covariance <- cov.fun.4(matrix(grid_p[Qmatrix_p[1,k]:Qmatrix_p[2,k],], qq[k], ncol(grid_p)),
                                         E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat)+ error                        
               }) 
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()
  Sigma_yij <- list()
  for(j in 1:(M-1)){
    N <- length(Neighbor_p[[j]])
    A1=A <- Neighbor_p[[j]]
    for(i in 1:N){
      A <- unique(c(A, Neighbor_p[[A1[i]]])) 
    }
    A <- sort(A)   
    Sigma_yij[[j]] <- lapply((j+1):M, function(i){ 
                        if(sum(i == A)){ 
                          covariance <- cov.fun.3(matrix(grid_p[Qmatrix_p[1,j]:Qmatrix_p[2,j],],(qq[j]),ncol(grid_p)),
                                                  matrix(grid_p[Qmatrix_p[1,i]:Qmatrix_p[2,i],],(qq[i]),ncol(grid_p)),
                                                  E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat)
                        }else{
                          return(0)
                        }                      
                      })
  }
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()  
  Sigma_yii_inverse <- lapply(1:M, function(i){
                         chol2inv( chol(Sigma_yii[[i]]) ) 
                       })   
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()
  predictsize <- nrow(grids)
  grids_p_region <- region.fun(grids, center_p)
  grids_label <- sapply(1:predictsize, function(i){ sum( grids_p_region[1:i]==grids_p_region[i] ) })
  pcov <- list()
  for(u in 1:M){
      #print(u)
      o <- which(grids_p_region==u)
      if(length(o)==0){
        pcov[[u]] <- 0
      }else{
        A <- sort(c(u,Neighbor_p[[u]]))
        L <- length(Neighbor_p[[u]])+1   
        pcov[[u]] <- lapply(1:L, function(i){
                       cov.fun.3( matrix(grids[o,],length(o),ncol(grids)), 
                                  matrix(grid_p[Qmatrix_p[1,A[i]]:Qmatrix_p[2,A[i]],],qq[A[i]],ncol(grids)), 
                                  E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat)
                     })                      
      }   
  }           
  ttt <- proc.time()-tt 
  #print(ttt)

  tt <- proc.time()  
  R_pred <- sapply(1:predictsize, function(k){
    u <- grids_p_region[k]
    A <- sort(c(u,Neighbor_p[[u]]))
    L <- length(Neighbor_p[[u]])+1                   
    bi <- lapply(1:L, function (i) { matrix(pcov[[ u ]][[ i ]][grids_label[k],] %*% Sigma_yii_inverse[[A[i]]], qq[A[i]], 1) } )
    hat_z_local <- matrix( sapply(1:L, function (i) { crossprod(bi[[i]] , matrix(Z_p[Qmatrix_p[1,A[i]]:Qmatrix_p[2,A[i]]], qq[A[i]], 1)) }), L, 1)
    di <- matrix( sapply(1:L, function (i) { tcrossprod( t(bi[[i]]) , pcov[[ u ]][[ i ]][grids_label[k],] ) }), 1, L )
    MM <- matrix(, L, L)
    for (i in 1:L) {
      for (j in 1:L) {
        if(i==j){ 
           MM[i,j] <- crossprod( bi[[i]] , Sigma_yii[[A[i]]] %*% bi[[j]] )
        }else if (A[i]>A[j]){
           MM[i,j] <- crossprod( bi[[i]] , t(Sigma_yij[[A[j]]][[A[i]-A[j]]]) %*% bi[[j]] )
        }else if (A[i]<A[j]){
           MM[i,j] <- crossprod( bi[[i]] , Sigma_yij[[A[i]]][[A[j]-A[i]]] %*% bi[[j]] )
        }  
      }
    }
    pred <- di %*% chol2inv(chol(MM)) %*% hat_z_local
    #va <- cov.fun.4(matrix(grids[k,],1,ncol(grids)),
    #                E$center, E$border, E$sigmasqvector_hat, E$alphavector_hat, E$nuvector_hat, E$a_hat) 
    #pred_sd <- sqrt(va - di %*% chol2inv(chol(MM)) %*% t(di))
    #return( c(pred, pred_sd) )
  })
  ttt <- proc.time()-tt
  #print(ttt)
  return( list(grids = grids, predictor = matrix(R_pred,predictsize,1) ))
}






STADK <- function(E, grid = E$grid, Z = E$Z, center_p, grids){
  
  Z_p <- order.fun(Z, grid, center_p)                  
  grid_p <- order.fun(grid, grid, center_p) 
  grid_p_region <- region.fun(grid_p, center_p)
  M <- nrow(center_p)                       
  Qmatrix_p <- matrix(,2,M)
  qq <- sapply(1:M, function(k){sum(grid_p_region==k)})
  Qmatrix_p[1,] <- cumsum(qq)-qq+1
  Qmatrix_p[2,] <- cumsum(qq)

  tt <- proc.time()   
  Sigma_yii <- lapply(1:M, function(k){
                 d <- dist1( matrix(grid_p[Qmatrix_p[1,k]:Qmatrix_p[2,k],], qq[k], ncol(grid_p)) )
                 error <- diag(E$tau_hat, qq[k])
                 covariance <- E$sigmasq_hat*matern(d, 1/E$alpha_hat, 0.5) + error                        
               }) 
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()
  Sigma_yij <- list()
  for(j in 1:(M-1)){   
    Sigma_yij[[j]] <- lapply((j+1):M, function(i){
                          d <- dist2( matrix(grid_p[Qmatrix_p[1,j]:Qmatrix_p[2,j],],(qq[j]),ncol(grid_p)),
                                      matrix(grid_p[Qmatrix_p[1,i]:Qmatrix_p[2,i],],(qq[i]),ncol(grid_p)) )
                          covariance <- E$sigmasq_hat*matern(d, 1/E$alpha_hat, 0.5) }) 
  }
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()   
  Sigma_yii_inverse <- lapply(1:M, function(i){
                         chol2inv( chol(Sigma_yii[[i]]) ) 
                       })   
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()
  predictsize <- nrow(grids)
  grids_p_region <- region.fun(grids, center_p)
  grids_label <- sapply(1:predictsize, function(i){ sum( grids_p_region[1:i]==grids_p_region[i] ) })
  pcov <- list()
  for(u in 1:M){
      #print(u)
      o <- which(grids_p_region==u)
      if(length(o)==0){
        pcov[[u]] <- 0
      }else{   
        pcov[[u]] <- lapply(1:M, function(i){
                       d <- dist2( matrix(grids[o,],length(o),ncol(grids)), 
                                   matrix(grid_p[Qmatrix_p[1,i]:Qmatrix_p[2,i],],qq[i],ncol(grids)) )
                       covariance <- E$sigmasq_hat*matern(d, 1/E$alpha_hat, 0.5)
                     })                      
      }   
  }           
  ttt <- proc.time()-tt 
  #print(ttt)
  
  tt <- proc.time()  
  R_pred <- sapply(1:predictsize, function(k){
    u <- grids_p_region[k]
    bi <- lapply(1:M, function (i) { matrix(pcov[[ u ]][[ i ]][grids_label[k],] %*% Sigma_yii_inverse[[ i ]], qq[i], 1) } )
    hat_z_local <- matrix( sapply(1:M, function (i) { crossprod(bi[[i]] , matrix(Z_p[Qmatrix_p[1,i]:Qmatrix_p[2,i]], qq[i], 1)) }), M, 1)
    di <- matrix( sapply(1:M, function (i) { tcrossprod( t(bi[[i]]) , pcov[[ u ]][[ i ]][grids_label[k],] ) }), 1, M )
    MM <- matrix(, M, M)
    for (i in 1:M) {
      for (j in 1:M) {
        if(i==j){ 
           MM[i,j] <- crossprod( bi[[i]] , Sigma_yii[[i]] %*% bi[[j]] )
        }else if (i>j){
           MM[i,j] <- crossprod( bi[[i]] , t(Sigma_yij[[j]][[i-j]]) %*% bi[[j]] )
        }else if (i<j){
           MM[i,j] <- crossprod( bi[[i]] , Sigma_yij[[i]][[j-i]] %*% bi[[j]] )
        }  
      }
    }

    pred <- di %*% chol2inv(chol(MM)) %*% hat_z_local
    #va <- E$sigmasq_hat*matern(0, 1/E$alpha_hat, 0.5) 
    #pred_sd <- sqrt(va - di %*% chol2inv(chol(MM)) %*% t(di))
    #return( c(pred, pred_sd) )
  })                       
  ttt <- proc.time()-tt
  #print(ttt)
  return( list(grids = grids, predictor = matrix(R_pred,predictsize,1) ))
}




STALDK <- function(E, grid = E$grid, Z = E$Z, center_p, grids){

  Z_p <- order.fun(Z, grid, center_p)                  
  grid_p <- order.fun(grid, grid, center_p) 
  grid_p_region <- region.fun(grid_p, center_p)
  M <- nrow(center_p)                       
  Qmatrix_p <- matrix(,2,M)
  qq <- sapply(1:M, function(k){sum(grid_p_region==k)})
  Qmatrix_p[1,] <- cumsum(qq)-qq+1
  Qmatrix_p[2,] <- cumsum(qq)

  vor_dxy_p <- deldir(center_p[,1], center_p[,2], rw = E$border)
  neighbor_p <- neighbor.fun(center_p, E$border)
  Neighbor_p <- lapply(1:M, function(i){
    A <- which( neighbor_p[,i]!=0 )
    B <- rbind(as.matrix(vor_dxy_p$dirsgs[A,1:2]),as.matrix(vor_dxy_p$dirsgs[A,3:4]))
    D <- dist1(B)
    del <- c()
    for(j in 1:dim(B)[1]){
      same <- which(D[j,]==0)
      del <-c( del,same[which((same-j)>0)] )
    }
    if(length(del)!=0){
        B <- B[-del, ]
    }
    E1 <- dist2(B,as.matrix(vor_dxy_p$dirsgs[,1:2]))
    E2 <- dist2(B,as.matrix(vor_dxy_p$dirsgs[,3:4]))
    F1 <- sapply(1:dim(neighbor_p)[1], function(j){prod(E1[,j])==0})
    F2 <- sapply(1:dim(neighbor_p)[1], function(j){prod(E2[,j])==0})
    relation <- which(F1+F2 != 0)
    G <- sapply(1:M, function(k){sum(neighbor_p[relation,k])})
    return( setdiff(which( G!=0 ),i) )
  })

  tt <- proc.time()  
  Sigma_yii <- lapply(1:M, function(k){
                 d <- dist1( matrix(grid_p[Qmatrix_p[1,k]:Qmatrix_p[2,k],], qq[k], ncol(grid_p)) )
                 error <- diag(E$tau_hat, qq[k])
                 covariance <- E$sigmasq_hat*matern(d, 1/E$alpha_hat, 0.5) + error                            
               }) 
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()
  Sigma_yij <- list()
  for(j in 1:(M-1)){
    N <- length(Neighbor_p[[j]])
    A1=A <- Neighbor_p[[j]]
    for(i in 1:N){
      A <- unique(c(A, Neighbor_p[[A1[i]]])) 
    }
    A <- sort(A) 
    Sigma_yij[[j]] <- lapply((j+1):M, function(i){ 
                        if(sum(i == A)){
                          d <- dist2( matrix(grid_p[Qmatrix_p[1,j]:Qmatrix_p[2,j],],(qq[j]),ncol(grid_p)),
                                      matrix(grid_p[Qmatrix_p[1,i]:Qmatrix_p[2,i],],(qq[i]),ncol(grid_p)) )
                          covariance <- E$sigmasq_hat*matern(d, 1/E$alpha_hat, 0.5)
                        }else{
                          return(0)
                        }                      
                      })
  }
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()   
  Sigma_yii_inverse <- lapply(1:M, function(i){
                         chol2inv( chol(Sigma_yii[[i]]) ) 
                       })   
  ttt <- proc.time()-tt
  #print(ttt)

  tt <- proc.time()
  predictsize <- nrow(grids)
  grids_p_region <- region.fun(grids, center_p)
  grids_label <- sapply(1:predictsize, function(i){ sum( grids_p_region[1:i]==grids_p_region[i] ) })
  pcov <- list()
  for(u in 1:M){
      #print(u)
      o <- which(grids_p_region==u)
      if(length(o)==0){
        pcov[[u]] <- 0
      }else{
        A <- sort(c(u,Neighbor_p[[u]]))
        L <- length(Neighbor_p[[u]])+1
        pcov[[u]] <- lapply(1:L, function(i){
                       d <- dist2( matrix(grids[o,],length(o),ncol(grids)), 
                                   matrix(grid_p[Qmatrix_p[1,A[i]]:Qmatrix_p[2,A[i]],],qq[A[i]],ncol(grids)) )
                       covariance <- E$sigmasq_hat*matern(d, 1/E$alpha_hat, 0.5)
                     } )                     
      }   
  }           
  ttt <- proc.time()-tt 
  #print(ttt)

  tt <- proc.time()  
  R_pred <- sapply(1:predictsize, function(k){
    u <- grids_p_region[k]
    A <- sort(c(u,Neighbor_p[[u]]))
    L <- length(Neighbor_p[[u]])+1                   
    bi <- lapply(1:L, function (i) { matrix(pcov[[ u ]][[ i ]][grids_label[k],] %*% Sigma_yii_inverse[[A[i]]], qq[A[i]], 1) } )
    hat_z_local <- matrix( sapply(1:L, function (i) { crossprod(bi[[i]] , matrix(Z_p[Qmatrix_p[1,A[i]]:Qmatrix_p[2,A[i]]], qq[A[i]], 1)) }), L, 1)
    di <- matrix( sapply(1:L, function (i) { tcrossprod( t(bi[[i]]) , pcov[[ u ]][[ i ]][grids_label[k],] ) }), 1, L )
    MM <- matrix(, L, L)
    for (i in 1:L) {
      for (j in 1:L) {
        if(i==j){ 
           MM[i,j] <- crossprod( bi[[i]] , Sigma_yii[[A[i]]] %*% bi[[j]] )
        }else if (A[i]>A[j]){
           MM[i,j] <- crossprod( bi[[i]] , t(Sigma_yij[[A[j]]][[A[i]-A[j]]]) %*% bi[[j]] )
        }else if (A[i]<A[j]){
           MM[i,j] <- crossprod( bi[[i]] , Sigma_yij[[A[i]]][[A[j]-A[i]]] %*% bi[[j]] )
        }  
      }
    }
    pred <- di %*% chol2inv(chol(MM)) %*% hat_z_local
    #va <- E$sigmasq_hat*matern(0, 1/E$alpha_hat, 0.5)  
    #pred_sd <- sqrt(va - di %*% chol2inv(chol(MM)) %*% t(di))
    #return( c(pred, pred_sd) )
  })
  ttt <- proc.time()-tt
  #print(ttt)
  return( list(grids = grids, predictor = matrix(R_pred,predictsize,1) ))
}




