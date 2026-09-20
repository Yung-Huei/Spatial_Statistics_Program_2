
# setwd("C:/Users/Betty/Desktop/Table4")
load("DeepKriging_seed0/DeepKriging.RData")
MAE_test0 <- MAE_test
RMSE_test0 <- RMSE_test

load("DeepKriging_seed1/DeepKriging.RData")
MAE_test1 <- MAE_test
RMSE_test1 <- RMSE_test

load("DeepKriging_seed2/DeepKriging.RData")
MAE_test2 <- MAE_test
RMSE_test2 <- RMSE_test

load("DeepKriging_seed3/DeepKriging.RData")
MAE_test3 <- MAE_test
RMSE_test3 <- RMSE_test

load("DeepKriging_seed4/DeepKriging.RData")
MAE_test4 <- MAE_test
RMSE_test4 <- RMSE_test

# Five seeds
MAE_test <- c(
    MAE_test0,
    MAE_test1,
    MAE_test2,
    MAE_test3,
    MAE_test4
)

RMSE_test <- c(
    RMSE_test0,
    RMSE_test1,
    RMSE_test2,
    RMSE_test3,
    RMSE_test4
)

# Mean
MAE_test_ave <- mean(MAE_test)
RMSE_test_ave <- mean(RMSE_test)

# Sample SD
MAE_test_sd <- sd(MAE_test)
RMSE_test_sd <- sd(RMSE_test)

cat("\nTesting:\n")

cat("MAE  =", MAE_test_ave,
    " (SD =", MAE_test_sd, ")\n")

cat("RMSE =", RMSE_test_ave,
    " (SD =", RMSE_test_sd, ")\n")