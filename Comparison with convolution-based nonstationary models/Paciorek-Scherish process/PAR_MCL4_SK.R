# ============================================================
# 5.1.2 Comparison with convolution-based nonstationary modeling
#
# Paciorek & Schervish (2006)
#
# Data generating model:
#   True nonstationary covariance
#
# Fitted model:
#   Proposed nonstationary covariance
#   with estimated parameters
#
# Spatial partitioning:
#   SPATIV
#   Candidate numbers of partitions: k = 2, 3, 4, 5, 6
#   The number of partitions is selected by minimum BIC
#
# Parameter estimation:
#   MCL_exp()
#   Given the BIC-selected number of partitions
#
# Prediction:
#   Full simple kriging using estimated parameters
#
# Criteria:
#   MAE, MSPE, MSPE_Band, MSPE_Bandc, CRPS, INT, CVG
#   F_loss, KL_loss
#   Computational time
#
# Computational time:
#   1. Spatial partitioning + BIC selection
#   2. Parameter estimation
#   3. Simple kriging prediction
#
# Spatial grid: 60 x 60 = 3600 locations
# Observation sample: N = 1000
# Replications: T = 100
# ============================================================

library(MASS)
library(fields)

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

#library(convoSPAT)

load("code1e03new.RData")
load("spatdiv.RData")

# ============================================================
# 1. Simulation settings
# ============================================================

g <- 60
n <- g^2

T <- 100
N <- 1000

# True nugget variance
nugget <- 0.1

# True Paciorek-Schervish smoothness
smoothness <- 0.5

# Candidate numbers of spatial partitions
k_candidates <- c(2,3,4,5,6)

border = c(0,1,0,1)

# ============================================================
# 2. Spatial grid
# ============================================================

grid1 = grid2 = seq(
  1 / (2 * g),
  1 - 1 / (2 * g),
  length.out = g
)

G <- expand.grid(grid1, grid2)

grids <- matrix(
  cbind(G$Var1, G$Var2),
  nrow = g^2,
  ncol = 2
)


# ============================================================
# 3. Spatially varying range
#
# alpha(s) = 4 * 4^s1
# ============================================================

aRangeObj <- list()

class(aRangeObj) <- "PS_range"

predict.PS_range <- function(object, newdata, ...) {

  newdata <- as.matrix(newdata)

  4 * 4^newdata[, 1]
}


# ============================================================
# 4. Spatially varying variance
#
# sigma^2(s) = 0.25 * 4^s2
# ============================================================

rhoObj <- list()

class(rhoObj) <- "PS_variance"

predict.PS_variance <- function(object, newdata, ...) {

  newdata <- as.matrix(newdata)

  0.25 * 4^newdata[, 2]
}


# ============================================================
# 5. TRUE Paciorek-Schervish covariance
# ============================================================

C_obs <- Paciorek.cov(
  x1 = grids,
  x2 = grids,
  aRangeObj = aRangeObj,
  rhoObj = rhoObj,
  smoothness = smoothness
)

#dim(C_obs)

# ============================================================
# 6. Cholesky decomposition of TRUE covariance
# ============================================================

COV_0.5 <- t(chol(C_obs))

# ============================================================
# 7. Distance matrix for ALL 3600 grid locations
#
# Used for F_loss and KL_loss.
#
# d[i,j] = Euclidean distance between grid locations
# ============================================================

d <- as.matrix(
  dist(grids)
)


# ============================================================
# 8. Storage
# ============================================================

pred_SK <- matrix(
  NA,
  nrow = T,
  ncol = n
)

sd_SK <- matrix(
  NA,
  nrow = T,
  ncol = n
)

SPEu_SK <- matrix(
  NA,
  nrow = T,
  ncol = n
)

mspe_SK <- matrix(
  NA,
  nrow = T,
  ncol = 1
)

mae_SK <- matrix(
  NA,
  nrow = T,
  ncol = 1
)

crps_SK <- matrix(
  NA,
  nrow = T,
  ncol = 1
)

int_SK <- matrix(
  NA,
  nrow = T,
  ncol = 1
)

cvg_SK <- matrix(
  NA,
  nrow = T,
  ncol = 1
)


# ------------------------------------------------------------
# Covariance losses
# ------------------------------------------------------------

F_loss <- matrix(
  NA,
  nrow = T,
  ncol = 1
)

KL_loss <- matrix(
  NA,
  nrow = T,
  ncol = 1
)


# ------------------------------------------------------------
# Selected number of partitions
# ------------------------------------------------------------

k_hat_all <- numeric(T)


# ------------------------------------------------------------
# BIC for each candidate k
# ------------------------------------------------------------

BIC_all <- matrix(
  NA,
  nrow = T,
  ncol = length(k_candidates)
)

colnames(BIC_all) <- paste0(
  "k=",
  k_candidates
)


# ------------------------------------------------------------
# Computational time
#
# Column 1 = user
# Column 2 = system
# Column 3 = elapsed
#
# timePartition = spatial partitioning
# timeEstimation = MCL parameter estimation
# timeSK = prediction
# ------------------------------------------------------------

timePartition <- matrix(
  NA,
  nrow = T,
  ncol = 3
)

timeEstimation <- matrix(
  NA,
  nrow = T,
  ncol = 3
)

timeSK <- matrix(
  NA,
  nrow = T,
  ncol = 3
)

# ============================================================
# 9. Scoring functions
# ============================================================

# ------------------------------------------------------------
# CRPS
#
# x$mean = predictive mean
# x$sd   = predictive standard deviation
# y      = true value
# ------------------------------------------------------------

crps <- function(x, y) {

  z <- (y - x$mean) / x$sd

  scores <- x$sd * (
    z * (2 * pnorm(z) - 1) +
      2 * dnorm(z) -
      1 / sqrt(pi)
  )

  return(scores)
}


# ------------------------------------------------------------
# Interval Score
# alpha = 0.05 gives 95% prediction interval
# ------------------------------------------------------------

intscore <- function(x, y, alpha = 0.05) {

  hw <- -qnorm(alpha / 2) * x$sd

  scores <- 2 * hw +
    (2 / alpha) *
    (
      ((x$mean - hw) - y) *
        (y < x$mean - hw) +

      (y - (x$mean + hw)) *
        (y > x$mean + hw)
    )

  return(scores)
}


# ------------------------------------------------------------
# Coverage
# ------------------------------------------------------------

cvg <- function(x, y, alpha = 0.05) {

  hw <- -qnorm(alpha / 2) * x$sd

  scores <- (
    y >= x$mean - hw &
    y <= x$mean + hw
  )

  return(scores)
}

# ============================================================
# 10. Simulation
# ============================================================

for (u in 1:T) {

  cat("Replication:", u, "\n")

  set.seed(u - 1)

  # ==========================================================
  # 10.1 Generate TRUE latent spatial process
  #
  # y ~ N(0, COV)
  # ==========================================================

  y <- as.vector(
    COV_0.5 %*% rnorm(n, 0, 1)
  )


  # ==========================================================
  # 10.2 Add observation error
  #
  # z = y + epsilon
  #
  # epsilon ~ N(0, nugget)
  # ==========================================================

  z <- y +
    sqrt(nugget) *
    rnorm(n, 0, 1)


  # ==========================================================
  # 10.3 Select N = 1000 observations
  # ==========================================================

  sample400 <- sample(
    1:n,
    N,
    replace = FALSE
  )

  coords <- grids[
    sample400,
    ,
    drop = FALSE
  ]

  Z <- matrix(
    z[sample400],
    N,
    1
  )

  # ----------------------------------------------------------
  # Data frame for SPATIV
  # ----------------------------------------------------------

  coords2 <- data.frame(
    x = grids[sample400, 1],
    y = grids[sample400, 2]
  )

  # ==========================================================
  # 10.4 Spatial partitioning + BIC selection
  # ==========================================================

  BIC_k <- numeric(
    length(k_candidates)
  )

  center_list <- vector(
    "list",
    length(k_candidates)
  )


  # ----------------------------------------------------------
  # Start partitioning time
  # ----------------------------------------------------------

  t1 <- proc.time()

  for (jj in seq_along(k_candidates)) {

    kk <- k_candidates[jj]

    # ========================================================
    # 10.4.1 Obtain centers using spatdiv()
    # ========================================================

    
    PAR <- spatdiv(
        coor = coords2,
        data = Z,
        k = kk,
        area = 1
      )

    center_tmp <- PAR$center$coor
    center_tmp <- center_tmp[
      order(center_tmp[, 2]),
      ,
      drop = FALSE
    ]

    # --------------------------------------------------------
    # Store center corresponding to this k
    # --------------------------------------------------------

    center_list[[jj]] <- center_tmp

    BIC_k[jj] <- PAR$bic

  }

  # ----------------------------------------------------------
  # Partitioning time
  # ----------------------------------------------------------

  timePartition[u, ] <- (
    proc.time() - t1
  )[1:3]


  # ==========================================================
  # 10.5 Select k_hat based on minimum BIC
  # ==========================================================

  best <- which.min(
    BIC_k
  )

  k_hat <- k_candidates[best]

  k_hat_all[u] <- k_hat

  BIC_all[u, ] <- BIC_k


  # ----------------------------------------------------------
  # Select center corresponding to k_hat
  # ----------------------------------------------------------

  center <- as.matrix(center_list[[best]])

  # ==========================================================
  # 10.6 Parameter estimation using MCL_exp
  # ==========================================================

  t1 <- proc.time()

  ER <- MCL_exp(
    grid = coords,
    Z = Z,
    center = center,
    border = border
  )

  timeEstimation[u, ] <- (
    proc.time() - t1
  )[1:3]



  # ============================================================
  # 10.7 Estimated covariance
  # ============================================================

  t1 <- proc.time()

  B <- cov.fun.4(
    grids,
    center,
    border,
    ER$sigmasqvector_hat,
    ER$alphavector_hat,
    ER$nuvector_hat,
    ER$a_hat
  )



  # ============================================================
  # 10.11 Prediction covariance
  # ============================================================

  C_po <- B[
    ,
    sample400,
    drop = FALSE
  ]

  C_oo <- B[
    sample400,
    sample400,
    drop = FALSE
  ]


  # ============================================================
  # 10.12 Estimated observation covariance
  # ============================================================

  K_hat <- C_oo +
    ER$tau_hat * diag(N)


  # ============================================================
  # 10.13 Simple Kriging prediction
  # ============================================================

  # ------------------------------------------------------------
  # K^{-1} Z
  # ------------------------------------------------------------

  beta_hat <- solve(
    K_hat,
    Z
  )


  # ------------------------------------------------------------
  # Prediction mean
  # ------------------------------------------------------------

  pred <- as.vector(
    C_po %*% beta_hat
  )


  # ------------------------------------------------------------
  # Prediction variance
  # ------------------------------------------------------------

  Kinv_Cop <- solve(
    K_hat,
    t(C_po)
  )

  pred_var <- diag(B) -
    rowSums(
      C_po * t(Kinv_Cop)
    )


  # ------------------------------------------------------------
  # Numerical protection
  # ------------------------------------------------------------

  pred_var <- pmax(
    pred_var,
    1e-10
  )

  pred_sd <- sqrt(
    pred_var 
  )


  # ------------------------------------------------------------
  # Prediction time
  # ------------------------------------------------------------

  timeSK[u, ] <- (
    proc.time() - t1
  )[1:3]

  
  # ============================================================
  # 10.8 True covariance
  # ============================================================

  C <- C_obs


  # ============================================================
  # 10.9 Frobenius loss
  # ============================================================

  A <- B - C

  F_loss[u, ] <- sqrt(
    tr(A %*% t(A))
  )


  # ============================================================
  # 10.10 KL loss
  # ============================================================

  KL_loss[u, ] <- 0.5 * (
    tr(
      solve(B) %*% C
    ) -
      n +
      determinant(B)$modulus -
      determinant(C)$modulus
  )


  # ============================================================
  # 10.14 Store prediction and SD
  # ============================================================

  pred_SK[u, ] <- pred

  sd_SK[u, ] <- pred_sd


  # ============================================================
  # 10.15 Prediction error
  # ============================================================

  error <- pred - y

  SPEu_SK[u, ] <- error^2


  # ============================================================
  # 10.16 MAE
  # ============================================================

  mae_SK[u, ] <- mean(
    abs(error)
  )


  # ============================================================
  # 10.17 MSPE
  # ============================================================

  mspe_SK[u, ] <- mean(
    error^2
  )


  # ============================================================
  # 10.18 CRPS
  # ============================================================

  crps_SK[u, ] <- mean(
    crps(
      list(
        mean = pred,
        sd = pred_sd
      ),
      y
    ),
    na.rm = TRUE
  )


  # ============================================================
  # 10.19 95% Interval Score
  # ============================================================

  int_SK[u, ] <- mean(
    intscore(
      list(
        mean = pred,
        sd = pred_sd
      ),
      y,
      alpha = 0.05
    ),
    na.rm = TRUE
  )


  # ============================================================
  # 10.20 95% Coverage
  # ============================================================

  cvg_SK[u, ] <- mean(
    cvg(
      list(
        mean = pred,
        sd = pred_sd
      ),
      y,
      alpha = 0.05
    ),
    na.rm = TRUE
  )


  # ==========================================================
  # 10.21 Save
  # ==========================================================

  save.image(
    "PAR_MCL4_SK.RData"
  )
}

# ============================================================
# 11. Summary
# ============================================================

mean_MSPE <- mean(
  mspe_SK,
  na.rm = TRUE
)

se_MSPE <- sd(
  mspe_SK,
  na.rm = TRUE
) / sqrt(T)


mean_MAE <- mean(
  mae_SK,
  na.rm = TRUE
)

se_MAE <- sd(
  mae_SK,
  na.rm = TRUE
) / sqrt(T)


mean_CRPS <- mean(
  crps_SK,
  na.rm = TRUE
)

se_CRPS <- sd(
  crps_SK,
  na.rm = TRUE
) / sqrt(T)


mean_INT <- mean(
  int_SK,
  na.rm = TRUE
)

se_INT <- sd(
  int_SK,
  na.rm = TRUE
) / sqrt(T)


mean_CVG <- mean(
  cvg_SK,
  na.rm = TRUE
)

se_CVG <- sd(
  cvg_SK,
  na.rm = TRUE
) / sqrt(T)


mean_Floss <- mean(
  F_loss,
  na.rm = TRUE
)

se_Floss <- sd(
  F_loss,
  na.rm = TRUE
) / sqrt(T)


mean_KLloss <- mean(
  KL_loss,
  na.rm = TRUE
)

se_KLloss <- sd(
  KL_loss,
  na.rm = TRUE
) / sqrt(T)


# ============================================================
# 12. Computational time summary
# ============================================================

mean_timePartition <- apply(
  timePartition,
  2,
  mean,
  na.rm = TRUE
)

se_timePartition <- apply(
  timePartition,
  2,
  sd,
  na.rm = TRUE
) / sqrt(T)


mean_timeEstimation <- apply(
  timeEstimation,
  2,
  mean,
  na.rm = TRUE
)

se_timeEstimation <- apply(
  timeEstimation,
  2,
  sd,
  na.rm = TRUE
) / sqrt(T)


mean_timeSK <- apply(
  timeSK,
  2,
  mean,
  na.rm = TRUE
)

se_timeSK <- apply(
  timeSK,
  2,
  sd,
  na.rm = TRUE
) / sqrt(T)

mean_k_hat <- mean(
  k_hat_all,
  na.rm = TRUE
)

table_k_hat <- table(
  k_hat_all
)


# ============================================================
# 13. Print results
# ============================================================

cat("\n")
cat("====================================================\n")
cat("Scenario I\n")
cat("SPATIV Partition + MCL_exp + Simple Kriging\n")
cat("Estimated Parameters\n")
cat("====================================================\n\n")


cat("Prediction performance\n")
cat("----------------------\n")

cat(
  "MAE   = ",
  round(mean_MAE, 4),
  " (SE = ",
  round(se_MAE, 4),
  ")\n",
  sep = ""
)

cat(
  "MSPE  = ",
  round(mean_MSPE, 4),
  " (SE = ",
  round(se_MSPE, 4),
  ")\n",
  sep = ""
)

cat(
  "CRPS  = ",
  round(mean_CRPS, 4),
  " (SE = ",
  round(se_CRPS, 4),
  ")\n",
  sep = ""
)

cat(
  "INT95 = ",
  round(mean_INT, 4),
  " (SE = ",
  round(se_INT, 4),
  ")\n",
  sep = ""
)

cat(
  "CVG95 = ",
  round(mean_CVG, 4),
  " (SE = ",
  round(se_CVG, 4),
  ")\n",
  sep = ""
)


cat("\n")
cat("Covariance loss\n")
cat("----------------------\n")

cat(
  "F-loss = ",
  round(mean_Floss, 8),
  " (SE = ",
  round(se_Floss, 8),
  ")\n",
  sep = ""
)

cat(
  "KL-loss = ",
  round(mean_KLloss, 8),
  " (SE = ",
  round(se_KLloss, 8),
  ")\n",
  sep = ""
)


cat("\n")
cat("====================================================\n")
cat("Computational time\n")
cat("====================================================\n")


# ------------------------------------------------------------
# Spatial partitioning + BIC selection
# ------------------------------------------------------------

cat("\nSpatial partitioning + BIC selection:\n")
cat("--------------------------------------\n")

cat(
  "User    = ",
  round(mean_timePartition[1], 4),
  " (SE = ",
  round(se_timePartition[1], 4),
  ")\n",
  sep = ""
)

cat(
  "System  = ",
  round(mean_timePartition[2], 4),
  " (SE = ",
  round(se_timePartition[2], 4),
  ")\n",
  sep = ""
)

cat(
  "Elapsed = ",
  round(mean_timePartition[3], 4),
  " (SE = ",
  round(se_timePartition[3], 4),
  ")\n",
  sep = ""
)


# ------------------------------------------------------------
# Parameter estimation
# ------------------------------------------------------------

cat("\nParameter estimation (MCL_exp):\n")
cat("--------------------------------\n")

cat(
  "User    = ",
  round(mean_timeEstimation[1], 4),
  " (SE = ",
  round(se_timeEstimation[1], 4),
  ")\n",
  sep = ""
)

cat(
  "System  = ",
  round(mean_timeEstimation[2], 4),
  " (SE = ",
  round(se_timeEstimation[2], 4),
  ")\n",
  sep = ""
)

cat(
  "Elapsed = ",
  round(mean_timeEstimation[3], 4),
  " (SE = ",
  round(se_timeEstimation[3], 4),
  ")\n",
  sep = ""
)


# ------------------------------------------------------------
# Simple kriging prediction
# ------------------------------------------------------------

cat("\nSimple kriging prediction:\n")
cat("--------------------------\n")

cat(
  "User    = ",
  round(mean_timeSK[1], 4),
  " (SE = ",
  round(se_timeSK[1], 4),
  ")\n",
  sep = ""
)

cat(
  "System  = ",
  round(mean_timeSK[2], 4),
  " (SE = ",
  round(se_timeSK[2], 4),
  ")\n",
  sep = ""
)

cat(
  "Elapsed = ",
  round(mean_timeSK[3], 4),
  " (SE = ",
  round(se_timeSK[3], 4),
  ")\n",
  sep = ""
)

cat("\n")
cat("Selected number of partitions\n")
cat("-----------------------------\n")

cat(
  "Mean k_hat = ",
  round(mean_k_hat, 4),
  "\n",
  sep = ""
)

cat("\nFrequency of selected k:\n")
print(table_k_hat)
