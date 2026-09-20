# ============================================================
# 5.1.2 Comparison with convolution-based nonstationary modeling
#
# Paciorek & Schervish (2006)
#
# Data generating model:
#   True nonstationary covariance
#
# Fitted model:
#   Stationary exponential covariance + nugget
#
# Prediction:
#   Full simple kriging using estimated parameters
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

# True Paciorek-Schervish smoothness
smoothness <- 0.5


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
# Estimated parameters
#
# paraML:
#   [1] sigma^2
#   [2] 1 / alpha
#   [3] nugget
# ------------------------------------------------------------

paraML <- matrix(
  NA,
  nrow = T,
  ncol = 3
)

colnames(paraML) <- c(
  "sigma2",
  "inv_range",
  "nugget"
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

timeML <- matrix(
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
  # y ~ N(0, C_PS)
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
  # 10.4 Estimate stationary covariance parameters
  #
  # Model:
  #
  # C(h) = sigma^2 exp(-h / alpha)
  #
  # Parameterization:
  #
  # par[1] = sigma^2
  # par[2] = 1 / alpha
  # par[3] = nugget
  #
  # NOTE:
  # loglik.GRF2() must be available in your environment.
  # ==========================================================

  t1 <- proc.time()

  like <- function(par) {

    loglik.GRF2(
      coords = coords,
      data = Z,
      obj.model = NULL,
      cov.model = "exp",
      cov.pars = c(
        par[1],
        1 / par[2]
      ),
      nugget = par[3],
      kappa = 0.5,
      lambda = 1,
      psiR = 1,
      psiA = 0,
      trend = "cte",
      method.lik = "ML",
      compute.dists = TRUE,
      realisations = NULL
    )
  }


  E <- optim(
    par = c(
      var(Z),
      1,
      nugget
    ),
    fn = like,
    lower = c(
      1e-6,
      1e-6,
      1e-6
    ),
    upper = c(
      Inf,
      Inf,
      Inf
    ),
    method = "L-BFGS-B",
    control = list(
      fnscale = -1
    )
  )


  timeML[u, ] <- (
    proc.time() - t1
  )[1:3]


  paraML[u, ] <- E$par


  # ==========================================================
  # 10.5 Estimated stationary covariance over ALL 3600 points
  #
  # B = estimated covariance matrix
  #
  # B[i,j] =
  #   sigma_hat^2 *
  #   exp(-d_ij / alpha_hat)
  #
  # ==========================================================

  sigma2_hat <- E$par[1]

  inv_alpha_hat <- E$par[2]

  nugget_hat <- E$par[3]


  B <- sigma2_hat *
    exp(
      -inv_alpha_hat * d
    )


  # ==========================================================
  # 10.6 Covariance losses
  # ==========================================================

  A <- B - C_obs


  # Frobenius loss
  #
  # ||B - C||_F
  #

  F_loss[u, 1] <- sqrt(
    sum(A^2)
  )


  # ----------------------------------------------------------
  # KL divergence
  #
  # 0.5[
  #   tr(B^-1 C)
  #   - n
  #   + log|B|
  #   - log|C|
  # ]
  #
  # Use solve() rather than explicitly calculating inverse.
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


  # ==========================================================
  # 10.7 Estimated covariance at observed locations
  # ==========================================================

  d_oo <- as.matrix(
    dist(coords)
  )


  C_oo_hat <- sigma2_hat *
    exp(
      -inv_alpha_hat * d_oo
    )


  # ==========================================================
  # 10.8 Estimated covariance:
  #
  # prediction locations ALL 3600
  # versus observations
  # ==========================================================

  d_po <- fields::rdist(
    grids,
    coords
  )

  C_po_hat <- sigma2_hat *
    exp(
      -inv_alpha_hat * d_po
    )


  # ==========================================================
  # 10.9 Estimated covariance matrix for observations
  # + estimated nugget
  # ==========================================================

  K_hat <- C_oo_hat +
    nugget_hat * diag(N)


  # ==========================================================
  # 10.10 Estimated simple kriging
  #
  # beta = K_hat^{-1} Z
  # ==========================================================

  t1 <- proc.time()


  beta_hat <- solve(
    K_hat,
    Z
  )


  # ==========================================================
  # 10.11 Prediction mean
  # ==========================================================

  pred <- as.vector(
    C_po_hat %*% beta_hat
  )


  # ==========================================================
  # 10.12 Prediction variance
  #
  # Var[Y(s0) | Z]
  #
  # = C(s0,s0)
  #   - C_po K^-1 C_op
  #
  # Since we predict the latent spatial process y,
  # the diagonal of the covariance is sigma^2_hat.
  # ==========================================================

  Kinv_Cop <- solve(
    K_hat,
    t(C_po_hat)
  )


  pred_var <- sigma2_hat -
    rowSums(
      C_po_hat * t(Kinv_Cop)
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
  # 10.13 Computational time for simple kriging
  # ==========================================================

  timeSK[u, ] <- (
    proc.time() - t1
  )[1:3]


  # ==========================================================
  # 10.14 Store prediction and SD
  # ==========================================================

  pred_SK[u, ] <- pred

  sd_SK[u, ] <- pred_sd


  # ==========================================================
  # 10.15 Prediction error
  # ==========================================================

  error <- pred - y

  SPEu_SK[u, ] <- error^2


  # ==========================================================
  # 10.16 MAE
  # ==========================================================

  mae_SK[u, 1] <- mean(
    abs(error)
  )


  # ==========================================================
  # 10.18 MSPE
  # ==========================================================

  mspe_SK[u, 1] <- mean(
    error^2
  )


  # ==========================================================
  # 10.19 CRPS
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
  # 10.20 95% Interval Score
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
  # 10.21 95% Coverage
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
  # 10.22 Plot squared prediction error for u = 1
  # ==========================================================

  #if (u == 1) {

  #  SPEu_obj <- matrix(
  #    SPEu_SK[u, ],
  #    nrow = g,
  #    ncol = g,
  #    byrow = FALSE
  #  )

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
  # 10.23 Save
  # ==========================================================

  save.image(
    "MLE_SK.RData"
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

mean_timeML <- apply(
  timeML,
  2,
  mean,
  na.rm = TRUE
)

se_timeML <- apply(
  timeML,
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
# 13. Estimated parameter summary
# ============================================================

mean_paraML <- apply(
  paraML,
  2,
  mean,
  na.rm = TRUE
)

se_paraML <- apply(
  paraML,
  2,
  sd,
  na.rm = TRUE
) / sqrt(T)


# ============================================================
# 14. Print results
# ============================================================

cat("\n")
cat("====================================================\n")
cat("Paciorek-Schervish (2006)\n")
cat("Stationary Exponential Model + Simple Kriging\n")
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
  round(mean_Floss, 4),
  " (SE = ",
  round(se_Floss, 4),
  ")\n",
  sep = ""
)

cat(
  "KL-loss = ",
  round(mean_KLloss, 4),
  " (SE = ",
  round(se_KLloss, 4),
  ")\n",
  sep = ""
)


cat("\n")
cat("Estimated parameters\n")
cat("----------------------\n")

cat(
  "sigma^2 = ",
  round(mean_paraML[1], 4),
  " (SE = ",
  round(se_paraML[1], 4),
  ")\n",
  sep = ""
)

cat(
  "1/alpha = ",
  round(mean_paraML[2], 4),
  " (SE = ",
  round(se_paraML[2], 4),
  ")\n",
  sep = ""
)

cat(
  "alpha = ",
  round(1 / mean_paraML[2], 4),
  "\n",
  sep = ""
)

cat(
  "nugget = ",
  round(mean_paraML[3], 4),
  " (SE = ",
  round(se_paraML[3], 4),
  ")\n",
  sep = ""
)


cat("\n")
cat("Computational time\n")
cat("----------------------\n")

cat("ML estimation:\n")
print(round(mean_timeML, 4))

cat("\nSimple kriging:\n")
print(round(mean_timeSK, 4))

cat("\nSE of ML time:\n")
print(round(se_timeML, 4))

cat("\nSE of SK time:\n")
print(round(se_timeSK, 4))

