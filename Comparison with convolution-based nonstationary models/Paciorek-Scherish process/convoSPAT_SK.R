# ============================================================
# 5.1.2 Comparison with convolution-based nonstationary modeling
#
# Paciorek & Schervish (2006)
#
# Data generating model:
#   True nonstationary covariance
#
# Fitted model:
#   Convolution-based nonstationary spatial covariance
#   fitted using convoSPAT::NSconvo_fit()
#
# Covariance specification:
#   - exponential correlation
#   - locally isotropic spatially varying range
#   - spatially varying process variance
#   - constant nugget
#   - 16 mixture-component locations
#
# Parameter estimation:
#   Maximum likelihood (ML)
#
# Prediction:
#   Simple kriging using the fitted nonstationary covariance model
#
# Prediction evaluation:
#   MAE, MSPE, CRPS, Interval Score (INT), Coverage (CVG)
#
# Covariance evaluation:
#   Frobenius loss (F-loss)
#   Kullback-Leibler divergence (KL-loss)
#
# Computational evaluation:
#   ML fitting time
#   Simple kriging prediction time
#
# Spatial domain:
#   [0,1] x [0,1]
#
# Prediction grid:
#   60 x 60 = 3600 locations
#
# Observation sample:
#   N = 1000 randomly selected locations
#
# Number of replications:
#   T = 100
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

library(convoSPAT)

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

COV <- Paciorek.cov(
  x1 = grids,
  x2 = grids,
  aRangeObj = aRangeObj,
  rhoObj = rhoObj,
  smoothness = smoothness
)


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

#d <- as.matrix(
#  dist(grids)
#)

# ============================================================
# 8. Storage
# ============================================================

# ============================================================
# Estimated parameters from convoSPAT
# ============================================================

# Estimated intercept
intercept_SK <- matrix(
  NA,
  nrow = T,
  ncol = 1
)

# Estimated measurement error variance (nugget variance)
tausq_SK <- matrix(
  NA,
  nrow = T,
  ncol = 1
)

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
  # 10.4 Fit convolution-based nonstationary model
  #
  # convoSPAT specification:
  #   - 16 mixture-component centers
  #   - regular 4 x 4 grid over [0,1]^2
  #   - exponential correlation
  #   - locally isotropic
  #   - spatially varying process variance
  #   - spatially varying range
  #   - constant nugget
  #   - kappa = 0.5 (exponential correlation)
  # ==========================================================

  t1 <- proc.time()

  # ----------------------------------------------------------
  # Mixture-component centers
  # ----------------------------------------------------------

  mc.locations <- expand.grid(
    seq(0.125, 0.875, length.out = 4),
    seq(0.125, 0.875, length.out = 4)
  )

  mc.locations <- as.matrix(mc.locations)

  # ----------------------------------------------------------
  # Determine fit radius corresponding to approximately
  # 200 nearest observations
  # ----------------------------------------------------------

  dist_mc <- fields::rdist(
    mc.locations,
    coords
  )

  # Distance to the 200th nearest observation
  radius_200 <- apply(
    dist_mc,
    1,
    function(x) sort(x, partial = 200)[200]
  )

  # Use the largest radius so that every mixture component
  # has at least 200 observations available
  fit.radius <- max(radius_200)

  # Check local sample sizes
  mc.N <- convoSPAT::mc_N(
    coords = coords,
    mc.locations = mc.locations,
    fit.radius = fit.radius
  )

  #cat(
  #  "Local sample sizes:",
  #  paste(mc.N, collapse = ", "),
  #  "\n"
  #)

  #cat(
  #  "Minimum local sample size:",
  #  min(mc.N),
  #  "\n"
  #)

  # ----------------------------------------------------------
  # convoSPAT nonstationary model
  # ----------------------------------------------------------

  NSfit <- convoSPAT::NSconvo_fit(
  
    coords = coords,
    data = Z,
  
    cov.model = "exponential",
  
    mean.model = Z ~ 1,
  
    mc.locations = mc.locations,
  
    fit.radius = fit.radius,
  
    # spatially varying process variance
    ns.variance = TRUE,
  
    # constant nugget
    ns.nugget = FALSE,
  
    # locally isotropic covariance
    local.aniso = FALSE,
  
    # exponential correlation
    # nu = 0.5
    fix.kappa = TRUE,
    kappa = 0.5,
  
    # maximum likelihood
    method = "ml",
  
    print.progress = FALSE
  )

  # ----------------------------------------------------------
  # convoSPAT computational time
  # ----------------------------------------------------------

  timeML[u, ] <- (
    proc.time() - t1
  )[1:3]

  # ----------------------------------------------------------
  # Store estimated intercept and measurement error
  # ----------------------------------------------------------

  intercept_SK[u, 1] <- NSfit$beta.GLS[1]

  # tausq.est is the estimated nugget variance
  tausq_SK[u, 1] <- NSfit$tausq.est

  # ==========================================================
  # 10.5 Nonstationary kriging prediction
  # ==========================================================

  t1 <- proc.time()

  pred.NS <- predict(
    NSfit,
    pred.coords = grids
  )

  pred <- as.vector(
    pred.NS$pred.means
  )

  pred_sd <- as.vector(
    pred.NS$pred.SDs
  )

  # ----------------------------------------------------------
  # Remove measurement error variance
  #
  # pred.NS$pred.SDs includes the nugget variance.
  # The target y is the latent spatial process,
  # so CRPS / INT / CVG should use the predictive SD
  # for the latent process.
  # ----------------------------------------------------------

  pred_sd_latent <- sqrt(
    pmax(
      pred_sd^2 - NSfit$tausq.est,
      0
    )
  )

  # ----------------------------------------------------------
  # Computational time
  # ----------------------------------------------------------

  timeSK[u, ] <- (
    proc.time() - t1
  )[1:3]


  # ==========================================================
  # 10.6 Store prediction and SD
  # ==========================================================

  pred_SK[u, ] <- pred

  sd_SK[u, ] <- pred_sd_latent

  # ==========================================================
  # 10.7 Estimated covariance B at ALL 3600 grid locations
  #
  # Fitted model:
  #   - exponential covariance
  #   - locally isotropic
  #   - spatially varying process variance
  #   - constant nugget
  #
  # Based on the covariance construction used in NSconvo_fit()
  # ==========================================================

  # ----------------------------------------------------------
  # 2. Distance from prediction grids to mixture centers
  # ----------------------------------------------------------

  dist_grid_mc <- fields::rdist(
    grids,
    NSfit$mc.locations
  )

  # dimension should be 3600 x 16

  # ----------------------------------------------------------
  # 3. Mixture-component weights
  #
  # Same formula as NSconvo_fit():
  #
  #   exp(-d^2 / (2 * lambda.w))
  # ----------------------------------------------------------

  weights_unnorm <- exp(
    -dist_grid_mc^2 /
      (2 * NSfit$lambda.w)
  )

  weights_grid <- weights_unnorm /
    rowSums(weights_unnorm)

  # ----------------------------------------------------------
  # 4. Interpolate local kernel parameters
  #
  # Because local.aniso = FALSE:
  #
  #   kernel =
  #
  #       [ lambda_k      0     ]
  #       [     0      lambda_k ]
  #
  # Therefore only kernel[1,1] is needed.
  # ----------------------------------------------------------

  mc.lambda <- NSfit$mc.kernels[1, 1, ]

  grid.lambda <- as.vector(
    weights_grid %*% mc.lambda
  )

  # ----------------------------------------------------------
  # 6. Construct pairwise average local range
  #
  # For locally isotropic covariance:
  #
  #   Sigma_i = lambda_i I
  #
  # and
  #
  #   Sigma_ij = (Sigma_i + Sigma_j) / 2
  #
  # Therefore:
  #
  #   Sigma_ij = ((lambda_i + lambda_j)/2) I
  # ----------------------------------------------------------

  lambda_ij <- outer(
    grid.lambda,
    grid.lambda,
    FUN = "+"
  ) / 2


  # ----------------------------------------------------------
  # 8. Scaling factor
  #
  # Following NSconvo_fit():
  #
  # Scale.mat =
  #   sqrt(sqrt(det(Sigma_i)))
  #   *
  #   sqrt(1 / det(Sigma_ij))
  #   *
  #   sqrt(sqrt(det(Sigma_j)))
  #
  # For Sigma_i = lambda_i I:
  #
  # det(Sigma_i) = lambda_i^2
  #
  # This simplifies to:
  #
  #   Scale_ij =
  #   sqrt(lambda_i * lambda_j) / lambda_ij
  # ----------------------------------------------------------

  Scale.mat <- outer(
   sqrt(grid.lambda),
   sqrt(grid.lambda)
  ) / lambda_ij


  # ----------------------------------------------------------
  # 9. Pairwise Euclidean distances
  # ----------------------------------------------------------

  Dist.euclidean <- fields::rdist(
    grids,
    grids
  )


  # ----------------------------------------------------------
  # 10. Nonstationary transformed distance
  #
  # Since:
  #
  #   Sigma_ij^{-1}
  #       = (1 / lambda_ij) I
  #
  # therefore:
  #
  #   D_ij =
  #   ||s_i - s_j|| / sqrt(lambda_ij)
  # ----------------------------------------------------------

  Dist.mat <- Dist.euclidean /
    sqrt(lambda_ij)


  # ----------------------------------------------------------
  # 11. Exponential correlation
  #
  # cov.model = "exponential"
  #
  # rho(D) = exp(-D)
  # ----------------------------------------------------------

  Unscl.corr <- exp(
    -Dist.mat
  )


  # ----------------------------------------------------------
  # 12. Nonstationary correlation
  # ----------------------------------------------------------

  NS.corr <- Scale.mat *
    Unscl.corr


  # ----------------------------------------------------------
  # 13. Add spatially varying process variance
  #
  # Cov_ij =
  #
  # sqrt(sigma_i^2)
  # *
  # NS.corr_ij
  # *
  # sqrt(sigma_j^2)
  #
  # ----------------------------------------------------------
  # ------------------------------------------------------------
  # Process variance at the 16 mixture-component locations
  # ------------------------------------------------------------

  mc.sigmasq <- NSfit$MLEs.save$sigmasq

  # ------------------------------------------------------------
  # Interpolate process variance to all 3600 prediction locations
  # ------------------------------------------------------------

  grid.sigmasq <- as.vector(
    weights_grid %*% mc.sigmasq 
  )

  # ------------------------------------------------------------
  # Construct covariance
  # ------------------------------------------------------------

  B <- outer(
    sqrt(grid.sigmasq),
    sqrt(grid.sigmasq)
  ) * NS.corr

  
  # ==========================================================
  # 10.8 Covariance losses
  # ==========================================================

  A <- B - COV

  # ----------------------------------------------------------
  # Frobenius loss
  # ----------------------------------------------------------

  F_loss[u, 1] <- sqrt(
    sum(A^2)
  )


  # ----------------------------------------------------------
  # KL divergence
  # ----------------------------------------------------------

  chol_B <- chol(B)

  chol_C <- chol(COV)


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
        COV
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
  # 10.9 Prediction error
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
        sd = pred_sd_latent
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
        sd = pred_sd_latent
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
        sd = pred_sd_latent
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
    "convoSPAT_SK.RData"
  )
}

# ============================================================
# 11. Summary
# ============================================================

mean_intercept <- mean(
  intercept_SK,
  na.rm = TRUE
)

se_intercept <- sd(
  intercept_SK,
  na.rm = TRUE
) / sqrt(T)

min_intercept <- min(
  intercept_SK,
  na.rm = TRUE
)

max_intercept <- max(
  intercept_SK,
  na.rm = TRUE
)

range_intercept <- max_intercept - min_intercept

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



# ============================================================
# 13. Print results
# ============================================================

cat("\n")
cat("====================================================\n")
cat("Scenario I\n")
cat("Proposed Nonstationary Covariance + Stationary Exponential ML\n")
cat("Simple Kriging Using Estimated Parameters\n")
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


cat("\nML:\n")
print(
  round(
    mean_timeML,
    4
  )
)

cat("\nSE of ML time:\n")
print(
  round(
    se_timeML,
    4
  )
)

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

cat("\n")
cat("Estimated intercept\n")
cat("----------------------\n")

cat(
  "Mean      = ",
  round(mean_intercept, 4),
  "\n",
  sep = ""
)

cat(
  "SE        = ",
  round(se_intercept, 4),
  "\n",
  sep = ""
)

cat(
  "Min       = ",
  round(min_intercept, 4),
  "\n",
  sep = ""
)

cat(
  "Max       = ",
  round(max_intercept, 4),
  "\n",
  sep = ""
)

cat(
  "Range     = ",
  round(range_intercept, 4),
  "\n",
  sep = ""
)