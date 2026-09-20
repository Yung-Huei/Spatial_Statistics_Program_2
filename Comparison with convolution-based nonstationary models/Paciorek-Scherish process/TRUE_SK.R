# ============================================================
# 5.1.2 Comparison with convolution-based nonstationary modeling
#
# Paciorek & Schervish (2006)
# True covariance + full simple kriging
#
# Spatial grid: 60 x 60 = 3600 locations
# Observation sample: N = 1000
# Replications: T = 100
# ============================================================

library(MASS)
library(fields)

# ============================================================
# 1. Simulation settings
# ============================================================

g <- 60
n <- g^2

T <- 100
N <- 1000

# Nugget variance
nugget <- 0.1

# Matérn smoothness
smoothness <- 0.5


# ============================================================
# 2. Spatial grid
#
# Grid cell centers:
#
#   1/(2g), 3/(2g), ..., 1 - 1/(2g)
#
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

# Check
#dim(grids)
# 3600 2


# ============================================================
# 3. Spatially varying range
#
# alpha(s) = 4 * 4^s1
#
# s1 = 0     --> alpha = 4
# s1 = 1     --> alpha = 16
# ============================================================

aRangeObj <- list()

class(aRangeObj) <- "PS_range"

predict.PS_range <- function(object, newdata, ...) {

  newdata <- as.matrix(newdata)

  # First spatial coordinate
  4 * 4^newdata[, 1]
}


# ============================================================
# 4. Spatially varying variance
#
# sigma^2(s) = 0.25 * 4^s2
#
# s2 = 0     --> sigma^2 = 0.25
# s2 = 1     --> sigma^2 = 1
# ============================================================

rhoObj <- list()

class(rhoObj) <- "PS_variance"

predict.PS_variance <- function(object, newdata, ...) {

  newdata <- as.matrix(newdata)

  # Second spatial coordinate
  0.25 * 4^newdata[, 2]
}


# ============================================================
# 5. True Paciorek-Schervish covariance
#
# C_obs[i,j] = C(s_i, s_j)
#
# This is the TRUE covariance used for:
#   1. data generation
#   2. simple kriging
# ============================================================

C_obs <- fields::Paciorek.cov(
  x1 = grids,
  x2 = grids,
  aRangeObj = aRangeObj,
  rhoObj = rhoObj,
  smoothness = smoothness
)

#dim(C_obs)
# 3600 3600


# ============================================================
# 6. Cholesky decomposition of TRUE covariance
#
# Z ~ N(0, C_obs)
#
# C_obs = L L'
# ============================================================

COV_0.5 <- t(chol(C_obs))


# ============================================================
# 7. Storage
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
# 8. Simulation + True Simple Kriging
# ============================================================

for (u in 1:T) {

  cat("Replication:", u, "\n")

  # ----------------------------------------------------------
  # Reproducible seed
  # ----------------------------------------------------------

  set.seed(u - 1)


  # ----------------------------------------------------------
  # Generate TRUE latent spatial process
  #
  # Y_true ~ N(0, C_obs)
  #
  # ----------------------------------------------------------

  y <- as.vector(
    COV_0.5 %*% rnorm(n, 0, 1)
  )


  # ----------------------------------------------------------
  # Add nugget
  #
  # Z_obs = Y_true + epsilon
  #
  # epsilon ~ N(0, nugget)
  # ----------------------------------------------------------

  z <- y +
    sqrt(nugget) *
    rnorm(n, 0, 1)


  # ----------------------------------------------------------
  # Randomly select N = 1000 observations
  # ----------------------------------------------------------

  sample400 <- sample(
    1:n,
    N,
    replace = FALSE
  )


  # ----------------------------------------------------------
  # Observation coordinates
  # ----------------------------------------------------------

  coords <- grids[sample400, , drop = FALSE]


  # ----------------------------------------------------------
  # Observed response
  # ----------------------------------------------------------

  Z <- z[sample400]

  # ==========================================================
  # 10.5 Covariance losses
  #
  # Since fitted covariance is exactly TRUE covariance:
  #
  # B = C_true = COV
  # ==========================================================

  B <- C_obs


  A <- B - C_obs


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

  chol_C <- chol(C_obs)

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
        C_obs
      )
    )
  )

  KL_loss[u, 1] <- 0.5 * (
    trace_BinvC -
      n +
      logdet_B -
      logdet_C
  )


  # ----------------------------------------------------------
  # True covariance components
  #
  # C_oo:
  # covariance among observed locations
  #
  # C_po:
  # covariance between all prediction locations
  # and observed locations
  # ----------------------------------------------------------

  C_oo <- C_obs[
    sample400,
    sample400,
    drop = FALSE
  ]

  C_po <- C_obs[
    ,
    sample400,
    drop = FALSE
  ]


  # ----------------------------------------------------------
  # Simple kriging
  #
  # Target:
  #
  #   y(s)
  #
  # Observation:
  #
  #   Z = y(s_obs) + epsilon
  #
  # Therefore:
  #
  # E[Y(s) | Z]
  #
  # = C_po
  #   [C_oo + nugget I]^{-1}
  #   Z
  #
  # Mean is assumed to be ZERO.
  # ----------------------------------------------------------

  t1 <- proc.time()


  # ----------------------------------------------------------
  # Solve:
  #
  # [C_oo + nugget I] beta = Z
  #
  # Instead of explicitly calculating inverse matrix.
  # ----------------------------------------------------------

  K_obs <- C_oo +
    nugget * diag(N)

  beta <- solve(
    K_obs,
    Z
  )


  # ----------------------------------------------------------
  # Predict all 3600 locations
  # ----------------------------------------------------------

  pred <- as.vector(
    C_po %*% beta
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
    K_obs,
    t(C_po)
  )


  pred_var <- diag(C_obs) -
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

  # ----------------------------------------------------------
  # Computational time
  # ----------------------------------------------------------

  timeSK[u, ] <- (
    proc.time() - t1
  )[1:3]


  # ----------------------------------------------------------
  # Store predictions
  # ----------------------------------------------------------

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
