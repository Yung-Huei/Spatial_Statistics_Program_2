# ============================================================
# 5.1.1 Sensitivity to the modeling partition
#
# Scenario I
#
# Data generating model:
#   Proposed nonstationary covariance
#
# Fitted model:
#   Proposed nonstationary covariance
#   with estimated parameters
#
# Spatial partitioning:
#   TRUE partition

# Parameter estimation:
#   MCL_exp()
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
#   1. Parameter estimation
#   2. Simple kriging prediction
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
#library(permute)
#library(Rfast)
#library(fields)
#library(vegan)
#colMin = Rfast::colMins

#library(convoSPAT)

load("code1e03new.RData")
#load("spatdiv.RData")

# ============================================================
# 1. Simulation settings
# ============================================================

g <- 60
n <- g^2

T <- 100
N <- 1000

# True nugget variance
nugget <- 0.1

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
# 3. Proposed nonostationary covariance
# ============================================================

center <- matrix(c(0.25,0.75,0.25,0.75,0.25,0.25,0.75,0.75),4,2)
border = c(0,1,0,1)
K <- dim(center)[1]                                                 #the number of partition
sigmasqvector <- matrix(c(1,1,0.25,0.25),K,1)
alphavector <- matrix(c(16,4,16,4),K,1)
nuvector <- matrix(c(1/2,1/2,1/2,1/2),K,1)
a <- 0.01 
                
COV <- cov.fun.4(grids, center, border, sigmasqvector, alphavector, nuvector, a)


# ============================================================
# 4. Cholesky decomposition of TRUE covariance
# ============================================================

COV_0.5 <- t(chol(COV))

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
# Define two spatial classes based on Band
# ============================================================

Band <- apply(
  grids,
  1,
  function(s) {
    min(
      abs(s[1] - 0.5),
      abs(s[2] - 0.5)
    ) <= 0.05
  }
)

# Band = TRUE
# Bandc = FALSE

n_B  <- sum(Band)
n_Bc <- sum(!Band)

cat("Number of locations in Band  :", n_B, "\n")
cat("Number of locations in Band^c:", n_Bc, "\n")

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

mspe_SK_Band <- numeric(T)
mspe_SK_Bandc <- numeric(T)

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
# Computational time
#
# Column 1 = user
# Column 2 = system
# Column 3 = elapsed
#
# timeEstimation = MCL parameter estimation
# timeSK = prediction
# ------------------------------------------------------------

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

  C <- COV


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
  # 10.17.1 MSPE for Band
  # ============================================================

  mspe_SK_Band[u] <- mean(
    error[Band]^2,
    na.rm = TRUE
  )


  # ============================================================
  # 10.17.2 MSPE for B^c
  # ============================================================

  mspe_SK_Bandc[u] <- mean(
    error[!Band]^2,
    na.rm = TRUE
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
    "TruePAR_MCL4_SK.RData"
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


# ============================================================
# Summary: MSPE for Band
# ============================================================

mean_MSPE_Band <- mean(
  mspe_SK_Band,
  na.rm = TRUE
)

se_MSPE_Band <- sd(
  mspe_SK_Band,
  na.rm = TRUE
) / sqrt(T)


# ============================================================
# Summary: MSPE for Band^c
# ============================================================

mean_MSPE_Bandc <- mean(
  mspe_SK_Bandc,
  na.rm = TRUE
)

se_MSPE_Bandc <- sd(
  mspe_SK_Bandc,
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


# ============================================================
# 13. Print results
# ============================================================

cat("\n")
cat("====================================================\n")
cat("Scenario I\n")
cat("True Partition + MCL_exp + Simple Kriging\n")
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
  "MSPE Band   = ",
  round(mean_MSPE_Band, 4),
  " (SE = ",
  round(se_MSPE_Band, 4),
  ")\n",
  sep = ""
)

cat(
  "MSPE Band^c = ",
  round(mean_MSPE_Bandc, 4),
  " (SE = ",
  round(se_MSPE_Bandc, 4),
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
