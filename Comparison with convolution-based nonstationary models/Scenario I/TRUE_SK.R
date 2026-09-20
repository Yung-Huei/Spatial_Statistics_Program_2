# ============================================================
# 5.1.2 Comparison with convolution-based nonstationary modeling
#
# Scenario I
#
# Data generating model:
#   Proposed nonstationary covariance
#
# Fitted model:
#   True covariance
#
# Prediction:
#   Full simple kriging using known parameters
#
# Criteria:
#   MAE, MSPE, CRPS, INT, CVG
#   F_loss, KL_loss
#   Computational time
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
# Computational time
#
# Column 1 = user
# Column 2 = system
# Column 3 = elapsed
# ------------------------------------------------------------

#timeML <- matrix(
#  NA,
#  nrow = T,
#  ncol = 3
#)

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

  Z <- z[sample400]


  # ==========================================================
  # 10.4 TRUE covariance
  #
  # Scenario I:
  # fitted covariance = TRUE covariance
  #
  # No parameter estimation is needed.
  # ==========================================================

  C_true <- COV

  tau_hat <- nugget


  # ==========================================================
  # 10.5 Covariance losses
  #
  # Since fitted covariance is exactly TRUE covariance:
  #
  # B = C_true = COV
  # ==========================================================

  B <- C_true


  A <- B - C_true


  # ----------------------------------------------------------
  # Frobenius loss
  # ----------------------------------------------------------

  F_loss[u, 1] <- sqrt(
    sum(A^2)
  )


  # ----------------------------------------------------------
  # KL divergence
  #
  # KL(C_true || B)
  #
  # Since B = C_true:
  #
  # KL = 0
  # ----------------------------------------------------------

  chol_B <- chol(B)

  chol_C <- chol(C_true)

  logdet_B <- 2 * sum(
    log(diag(chol_B))
  )

  logdet_C <- 2 * sum(
    log(diag(chol_C))
  )

  trace_BinvC <- sum(
    diag(
      solve(
        B,
        C_true
      )
    )
  )

  KL_loss[u, 1] <- 0.5 * (
    trace_BinvC -
      n +
      logdet_B -
      logdet_C
  )


  # ==========================================================
  # 10.6 TRUE covariance:
  #
  # prediction locations ALL 3600
  # versus observations
  # ==========================================================

  C_po <- C_true[
    ,
    sample400,
    drop = FALSE
  ]


  # ==========================================================
  # 10.7 TRUE covariance among observations
  # ==========================================================

  C_oo <- C_true[
    sample400,
    sample400,
    drop = FALSE
  ]


  # ==========================================================
  # 10.8 Observation covariance
  #
  # Z = Y + epsilon
  #
  # Var(Z) = C_oo + nugget I
  # ==========================================================

  K_true <- C_oo +
    tau_hat * diag(N)


  # ==========================================================
  # 10.9 Simple Kriging
  #
  # beta = K_true^{-1} Z
  # ==========================================================

  t1 <- proc.time()

  beta_true <- solve(
    K_true,
    Z
  )


  # ==========================================================
  # 10.10 Prediction mean
  #
  # E[Y(s0) | Z]
  #
  # = C_po K^{-1} Z
  # ==========================================================

  pred <- as.vector(
    C_po %*% beta_true
  )


  # ==========================================================
  # 10.11 Prediction variance
  #
  # Var[Y(s0) | Z]
  #
  # = C(s0,s0)
  #   - C_po K^{-1} C_op
  #
  # IMPORTANT:
  # We predict latent Y, NOT noisy Z.
  #
  # Therefore:
  # diagonal of C_true = latent variance
  # ==========================================================

  Kinv_Cop <- solve(
    K_true,
    t(C_po)
  )


  pred_var <- diag(C_true) -
    rowSums(
      C_po * t(Kinv_Cop)
    )


  # Numerical protection
  pred_var <- pmax(
    pred_var,
    1e-10
  )


  pred_sd <- sqrt(
    pred_var
  )


  # ==========================================================
  # 10.12 Computational time
  # ==========================================================

  timeSK[u, ] <- (
    proc.time() - t1
  )[1:3]


  # ==========================================================
  # 10.13 Store prediction and SD
  # ==========================================================

  pred_SK[u, ] <- pred

  sd_SK[u, ] <- pred_sd


  # ==========================================================
  # 10.14 Prediction error
  # ==========================================================

  error <- pred - y

  SPEu_SK[u, ] <- error^2


  # ==========================================================
  # 10.15 MAE
  # ==========================================================

  mae_SK[u, 1] <- mean(
    abs(error)
  )


  # ==========================================================
  # 10.16 MSPE
  # ==========================================================

  mspe_SK[u, 1] <- mean(
    error^2
  )


  # ==========================================================
  # 10.17 CRPS
  # ==========================================================

  crps_SK[u, 1] <- mean(
    crps(
      list(
        mean = pred,
        sd = pred_sd
      ),
      y
    ),
    na.rm = TRUE
  )


  # ==========================================================
  # 10.18 95% Interval Score
  # ==========================================================

  int_SK[u, 1] <- mean(
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


  # ==========================================================
  # 10.19 95% Coverage
  # ==========================================================

  cvg_SK[u, 1] <- mean(
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
  # 10.20 Plot squared prediction error for u = 1
  # ==========================================================

  # if (u == 1) {

  #   SPEu_obj <- matrix(
  #     SPEu_SK[u, ],
  #     nrow = g,
  #     ncol = g,
  #     byrow = FALSE
  #   )

   #  fields::image.plot(
   #    x = grid1,
   #    y = grid2,
   #    z = SPEu_obj,
   #    xlab = "s1",
   #    ylab = "s2",
   #    main = "Squared Prediction Error: u = 1"
   #  )
   #}


  # ==========================================================
  # 10.21 Save
  # ==========================================================

  save.image(
    "TRUE_SK.RData"
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
cat("TRUE Nonstationary Covariance + Simple Kriging\n")
cat("Known True Parameters\n")
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
cat("Computational time\n")
cat("----------------------\n")

cat("\nSimple kriging:\n")
print(
  round(
    mean_timeSK,
    4
  )
)

cat("\nSE of SK time:\n")
print(
  round(
    se_timeSK,
    4
  )
)
