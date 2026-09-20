# ============================================================
# Five-fit average
# ============================================================
setwd("C:/Users/Betty/Desktop/Table4")
load("DeepKriging_seed0/DeepKriging.RData")
pred_test0 <- pred_test

load("DeepKriging_seed1/DeepKriging.RData")
pred_test1 <- pred_test

load("DeepKriging_seed2/DeepKriging.RData")
pred_test2 <- pred_test

load("DeepKriging_seed3/DeepKriging.RData")
pred_test3 <- pred_test

load("DeepKriging_seed4/DeepKriging.RData")
pred_test4 <- pred_test


# Average prediction
pred_test_ave <- (
    pred_test0 +
    pred_test1 +
    pred_test2 +
    pred_test3 +
    pred_test4
) / 5


# Testing MAE
MAE_test <- mean(
    abs(
        y_test - pred_test_ave
    )
)


# Testing RMSE
RMSE_test <- sqrt(
    mean(
        (y_test - pred_test_ave)^2
    )
)


cat("\nTesting:\n")
cat("MAE  =", MAE_test, "\n")
cat("RMSE =", RMSE_test, "\n")
