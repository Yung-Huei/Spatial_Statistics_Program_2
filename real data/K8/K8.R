#sat.temp
#beta0 = mean

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

load("code1e03new.RData")
load("spatdiv.RData")
load("AllSatelliteTemps.RData")

D <- subset(all.sat.temps,!is.na(TrueTemp)) 
D2 <- subset(D,!is.na(MaskTemp)) #training
D3 <- subset(D,is.na(MaskTemp))  #testing
rm(D)
rm(all.sat.temps)
coords <- matrix(cbind(D2$Lon, D2$Lat), nrow(D2), 2)
intercept <- mean(D2$MaskTemp)
Temp <- matrix(D2$MaskTemp, nrow(D2), 1) - intercept
area <- ( max(D2$Lon) - min(D2$Lon) ) * ( max(D2$Lat) - min(D2$Lat) )
border <- c(min(D2$Lon), max(D2$Lon), min(D2$Lat), max(D2$Lat))

#################
## determine K ##
#################

t1 <- proc.time()

K_hat <- 8

set.seed(0)
n=NROW(Temp)
fold=sample(1:n,n)
xis=rep(NA,n)
centers=NULL
for(k in 1:30){
    id=fold[round(n/30*(k-1)+1):round(n/30*k)]
    xis[id] = cal_xi(coor=coords[id,],data=Temp[id,],nn=6,area=area)
    centers = rbind(centers,as.matrix(spatdiv(coor=coords[id,],data=Temp[id,],k=K_hat,area=area)$center$coor))
}
res <- vorocluster(coor=coords, data=xis, k=K_hat, cslot=centers, fig.show = F)
center <- res$center$coor
res$bic

proc.time()-t1
K_hat
center

rm(xis)
rm(res)
rm(centers)
save.image("Partition.RData")


####################################################
## MCL2 based on exponential correlation function ##
####################################################

set.seed(0)

t1 <- proc.time()

ER <- MCL0_exp.2(grid=coords, Z=Temp, center=center, border=border)

t3 <- proc.time()-t1
t3
ER$sigmasqvector_hat
ER$alphavector_hat
ER$nuvector_hat 
ER$a_hat 
ER$tau_hat 
ER$cl

save.image("Estimate.RData")

##########################
## LDK with different M ##
##########################

t1 <- proc.time()

grid <- matrix(cbind(D3$Lon, D3$Lat),nrow(D3),2)
M <- 324
Lx <- (border[2]-border[1])
Ly <- (border[4]-border[3])
sp1 = seq(border[1] + Lx/2/18, border[2] - Lx/2/18, l = 18)
sp2 = seq(border[3] + Ly/2/18, border[4] - Ly/2/18, l = 18)
Sp <- expand.grid(sp1, sp2)
center_p <- as.matrix(Sp)
                  
grid_p <- order.fun(coords, coords, center_p) 
grid_p_region <- region.fun(grid_p, center_p)                 
qq <- sapply(1:M, function(k){sum(grid_p_region==k)})
center_p <- center_p[-which(qq<200),]
M <- nrow(center_p)

LDK81 <- LDK(E=ER, center_p=center_p, grids=grid)

t4 <- proc.time()-t1
t4
save.image("Predict.RData")

#########################################
## criteria: MAE, RMSE, CRPS, INT, CVG ##
#########################################

crps <- function(predlist,trueobs) {
  z <- as.numeric((trueobs - predlist$mean) / predlist$sd)
  scores <- predlist$sd * (z *(2 * pnorm(z, 0, 1) - 1) +
                             2 * dnorm(z, 0, 1) - 1/sqrt(pi))
  return(scores)
}

intscore <- function(x, y, alpha=0.05) {
  hw <- -qnorm(alpha/2) * x$sd
  scores <- 2 * hw + (2/alpha) * (((x$mean - hw) - y) * (y < x$mean - hw) +
                                    (y - (x$mean + hw)) * (y > x$mean + hw))
  return(scores)
}

cvg <- function(x, y, alpha=0.05) {
  hw <- -qnorm(alpha/2) * x$sd
  scores <- y >= (x$mean - hw) & y <= (x$mean + hw)
  return(scores)
}

TT <- matrix(D3$TrueTemp,nrow(D3),1)
MAE_LDK81 <- mean( abs(LDK81$predictor + intercept - TT), na.rm=TRUE )
RMSE_LDK81 <- sqrt( mean( (LDK81$predictor + intercept - TT)^2, na.rm=TRUE ))
CRPS_LDK81 <- mean( crps( list(mean=(LDK81$predictor + intercept), sd=LDK81$sd), TT ), na.rm=TRUE )
INT_LDK81 <- mean( intscore( list(mean=(LDK81$predictor + intercept), sd=LDK81$sd), TT ), na.rm=TRUE )
CVG_LDK81 <- mean( cvg( list(mean=(LDK81$predictor+ intercept), sd=LDK81$sd), TT ), na.rm=TRUE )

MAE_LDK81 
RMSE_LDK81
CRPS_LDK81 
INT_LDK81 
CVG_LDK81


save.image("Predict.RData")
