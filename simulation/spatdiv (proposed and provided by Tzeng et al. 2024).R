install.packages("permute")
install.packages("Rfast")
install.packages("vegan")
library(permute)
library(Rfast)
library(fields)
library(vegan)
colMin = Rfast::colMins 
#library(sp)

cal_xi <- function(geodata,coor=geodata$coor,data=geodata$data,nn,area,minPts=5){
    dmat <- as.matrix(dist(coor))
	n <- NROW(coor)
    output <- numeric(n)  #x#output <- c()
    n <- nrow(coor)
    maxd <- (nn*area/n/pi)^0.5  #n*(pi*maxd^2/area)>=nn=(minPts+1)
    for(i in 1:n){
        y_inrange <- data[-i][dmat[-i,i]<=maxd]
        if(length(y_inrange)<minPts){# control the least point for estimating xi
            output[i] <- NA #x#output <- c(output,NA)
            next
        }else{
            y_center <- data[i]
            yis <- abs(y_inrange-y_center)**0.5
            dis <- dmat[-i,i][dmat[-i,i]<=maxd]
            xi <- (mean(yis)-0)/mean(dis)
            output[i] <- xi #x#output <- c(output,xi)
        }
    }
    matrix(output,ncol=1)
}
# Plots -------------------------------------------------------------------
Hovmoller_pt <- function(dat,x,y,z,colnum=4,text=NULL,textsize=4,discrete=FALSE,collim=NULL){
    if(discrete){
        dat[,z] <- as.factor(dat[,z])
        p <- ggplot(dat) + # take data
            geom_point(aes_string(x = x, y = y, col = z)) + # plot
            ylab(y) + # add y label
            xlab(x) + # add x label
            theme_bw() + # change theme
            theme(panel.background = element_rect(fill="white")) # background color
        #p <- p + scale_color_manual(values = rainbow(length(unique(dat[,z]))))
    }else{
        p <- ggplot(dat) + # take data
            geom_point(aes_string(x = x, y = y, col = z)) + # plot
            ylab(y) + # add y label
            xlab(x) + # add x label
            theme_bw() + # change theme
            theme(panel.background = element_rect(fill="white"))# background color
        if(is.null(collim)){
            p <- p + scale_color_gradientn(colours = rainbow(colnum))
        }else{
            p <- p + scale_color_gradientn(colours = rainbow(colnum),limits=collim)
        }
    }
    if(!is.null(text)){
        p <- p + annotate("text", x=Inf, y = Inf, col="black", vjust=1.2, hjust=1.2,
                          label=text,size=textsize)
    }
    p + theme(aspect.ratio=1)# force plot to be square shaped
}



# Clustering --------------------------------------------------------------

## deterministic k-means 
dkmCenter <- function(mydata, k){

  #################################################################################
  # computing the distance between every pair of points
  #################################################################################
  M <- as.matrix(dist(mydata))
  M <- M / max(M)   # min-max normalization of distances

  #################################################################################
  # MST heuristic for eps calculation
  #################################################################################
  mstedglengs <- spantree(dist(mydata))$dist
  eps <-  min(max(mstedglengs), 3.0 * IQR(mstedglengs) + quantile(mstedglengs, probs = c(0.75)))

  #################################################################################
  # rho computation:
  #################################################################################
  n <- nrow(mydata)
  rho <- vector("double", n)
  for(i in 1:n){
    epsnbrsi <- setdiff(which(M[ ,i] <= eps), i) #epsilon neghbors of data point i
    for(j in epsnbrsi){
      rho[i] <- rho[i] + exp(-1 * M[i, j] / eps)
    }
  }

  #################################################################################
  # min-max normalization of rho
  #################################################################################
  maxrho <- max(rho)
  minrho <- min(rho)
  rho <- rho - minrho
  rho <- rho / (maxrho - minrho)

  #################################################################################
  # Finding k initial centroids
  #################################################################################
  p <- which.max(rho) # highest density point
  prospect <- vector("double", n)
  C <- p                                        # C = {p}
  mind2cent <- rep(Inf, n)
  while(length(C) < k){
    for(j in 1:n){
      mind2cent[j] <- min(M[j, p], mind2cent[j])
      prospect[j] <- rho[j] * mind2cent[j]
    }
    p <- which.max(prospect) # point with maximum prospect = a centroid
    C <- union(C, p) # Added center to C
  }

  #################################################################################
  # call kmeans(Hartigan-Wong) with initial centroids in C
  #################################################################################
  #return(kmeans(mydata, mydata[C, ], iter.max = 50)$cluster)
  return(mydata[C, ])
}



## calculate distance between multiple points
distFun <- function(points1, points2) {# row is pt1    #x#
    points1 <- as.matrix(points1);points2 <- as.matrix(points2)
    fields::rdist(points1,points2)
}

## cluster data according to their distance to centers
get_cluster <- function(coor,center_coor){
    disttocenter <- distFun(as.matrix(coor),center_coor)
    cluster <- colMin(t(disttocenter*1.0))
    cluster
}
## calculate negative log-likilihood given centers
SD=function(x) sqrt(mean((x-mean(x))^2))
cal_n_loglik <- function(coor,data,center_coor){
    cluster <- get_cluster(coor,center_coor)
    k <- nrow(center_coor)# cluster number
    loglik <- 0
    for(i in 1:k){
        data_i <- data[cluster==i]
        if(length(data_i)>1){
            #loglik <- loglik + sum(log(dnorm(data_i,mean(data_i),SD(data_i))))
			ss=SD(data_i)
			loglik <- loglik + sum(log(dnorm((data_i-mean(data_i))/ss))-log(ss))
        }else{
            loglik <- loglik + -Inf
        }
    }
    return(-loglik)
}

# search cluster center #x#  main changes of cslot
search_center <- function(coor=coor,data=data,k, cslot=NULL){
    if(is.null(cslot)) cslot=coor
    n <- nrow(cslot)
    nloglik.now <- Inf    
	
	center_coor <- dkmCenter(cslot, k)  #x#
	nloglik.now <- cal_n_loglik(coor,data,center_coor) #x#	
	center <- list(coor=center_coor)  #x#
	
    # update center
    while(TRUE){
        cluster <- get_cluster(coor,center_coor)
        unchanged_cluster <- 0
        for(i in 1:k){
            nbind <- which(cluster==i)
            nb_coor <- cslot[nbind,]
            nloglik.next <- apply(nb_coor,1,function(x){cal_n_loglik(coor,data,rbind(center_coor[-i,],x))})
            if(min(nloglik.next)<nloglik.now){
                nloglik.now <- min(nloglik.next)
                newind=which.min(nloglik.next) #x#
                center_coor[i,] <- nb_coor[newind,]
                center$coor[i,] <- nb_coor[newind,]
                #x# center$data[i] <- data[nbind][newind]
                break
            }else{
                unchanged_cluster <- unchanged_cluster + 1
            }
        }
        if(unchanged_cluster==k) break
    }
    center
}

# calculate BIC
cal_bic <- function(coor,data,center_coor){
    if(is.list(center_coor) & (length(center_coor)==1))
      center_coor <- center_coor[[1]]    
    coor.rmNA <- coor[!is.na(data),]
    data.rmNA <- data[!is.na(data)]
    k <- nrow(center_coor)
    n <- nrow(coor.rmNA)
    cluster <- get_cluster(coor=coor.rmNA,center_coor = center_coor)
    nloglik <- cal_n_loglik(coor=coor.rmNA,data=data.rmNA,center_coor = center_coor)
    #x# bic <- 2*nloglik + 2*k*log(n)
	bic <- 2*nloglik + 4*k*log(n)
    bic
}

# spatial clustering
vorocluster <- function(coor, data ,k,cslot=NULL,fig.show=FALSE){
    coor.rmNA <- coor[!is.na(data),]
    data.rmNA <- data[!is.na(data)]
    if(nrow(coor.rmNA)<(2*k)) stop(paste("total",nrow(coor.rmNA) ,"points, there are not enough points for clustering"))
    center_res.rmNA <- search_center(coor=coor.rmNA,data=data.rmNA,k=k)
    cluster <- get_cluster(coor=coor,center_coor = center_res.rmNA$coor)
    if(fig.show){
        p <- Hovmoller_pt(cbind(coor,c=cluster),"x","y","c",discrete=TRUE) +
             geom_point(data=center_res.rmNA$coor,aes(x=x,y=y))
        print(p)
    }
    BIC=cal_bic(coor,data,center_res.rmNA)
    list(cluster=cluster,center=center_res.rmNA, bic=BIC)
}

spatdiv <- function(coor, data, k, area, nn=6, cslot=NULL, plot=FALSE) 
{
xi <- cal_xi(coor=coor,data=data,nn=nn,area=area)
result <- vorocluster(coor=coor,data=xi,k=k, cslot, fig.show = plot)
return(result)
}

##==========================================================
## example using simulation
##==========================================================

library(mvtnorm)
library(geoR)
library(ggplot2)
#clustering_sim <- function(n,lambdas,a,nn,fig.show=FALSE){
matern_cov <- function(distmat,sigmasq,nu,beta,tausq){
    # u = ||si-sj||
    # phi = beta = range parameter
    # kappa = nu = smoothness parameter
    covmat <- geoR::matern(u=distmat, phi=beta, kappa=nu)
    covmat <- covmat*sigmasq
    covmat[distmat==0] <- covmat[distmat==0] + tausq
    covmat
}    



n=500
a <- 0.25
k <- 2
lambdas <- c(0.1,0.5)
xrange <- c(0,1)
yrange <- c(0,1)
coor <- as.data.frame(cbind(x=runif(n,xrange[1],xrange[2]),y=runif(n,yrange[1],yrange[2])))
cluster <- rep(1,nrow(coor))
cluster[coor$x<=coor$y] <- 2
c1=matern_cov(as.matrix(dist(coor[cluster==1,])),3,nu=0.5,beta=lambdas[1],tausq=a)
c2=matern_cov(as.matrix(dist(coor[cluster==2,])),2,nu=1,beta=lambdas[2],tausq=a)
y=coor[,1]*NA
y[cluster==1]=c(mvtnorm::rmvnorm(1,sigma=c1))
y[cluster==2]=c( mvtnorm::rmvnorm(1,sigma=c2))
spatdiv(coor,data=y,k=2,area=1,plot=T)



# 主程式 spatdiv 說明
# input
#  coor:    二維座標 (nx2矩陣)
#  data:    對應 coor 位置上的觀測值 (長度n向量)
#  k:       要切成幾塊
#  area:    coor 所在domain的面積大小 
#  nn:      turning parameter 預設為 6 (透過周圍多少點判斷差異)
#  cslot:   限制搜尋各分塊的中心點集合 (mx2矩陣)，預設為NULL，即使用 coor 
#  plot:    是否要畫圖看一下分割結果，預設為 FALSE
#
# output
#  cluster: 根據 center 做 Voronoi Diagram 分割，coor每點對應的子區域
#  center:  每塊用來分割的中心點
#  bic:     分割後的 Bayesian Information Criterion  