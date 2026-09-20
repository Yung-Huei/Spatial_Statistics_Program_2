# ============================================================
# 2D Wendland DeepKriging
#
# CPU version
#
# Output for Python uncertainty estimation:
#
#   y_train.csv
#   coords_train.csv
#   pred_train.csv
#   coords_test.csv
#   pred_test.csv
#   y_test.csv
#
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

library(keras3)
library(geoR)


# ============================================================
# 1. Settings
# ============================================================

# -------------------------
# Random seed
# -------------------------

set.seed(0)


# -------------------------
# Number of CPU cores
# -------------------------

num_cores <- 49

cat("Requested CPU cores:", num_cores, "\n")


# ============================================================
# 2. Configure TensorFlow CPU threads
# ============================================================

Sys.setenv(
    OMP_NUM_THREADS = num_cores,
    TF_NUM_INTRAOP_THREADS = num_cores,
    TF_NUM_INTEROP_THREADS = num_cores
)


# ============================================================
# 3. Load data
# ============================================================

load("AllSatelliteTemps.RData")


D <- subset(
    all.sat.temps,
    !is.na(TrueTemp)
)


# Training data
D2 <- subset(
    D,
    !is.na(MaskTemp)
)


# Testing data
D3 <- subset(
    D,
    is.na(MaskTemp)
)


rm(D)
rm(all.sat.temps)


# ============================================================
# 4. Coordinates
# ============================================================

t1 <- proc.time()

coords_train <- as.matrix(
    cbind(
        D2$Lon,
        D2$Lat
    )
)

coords_test <- as.matrix(
    cbind(
        D3$Lon,
        D3$Lat
    )
)

colnames(coords_train) <- c(
    "Lon",
    "Lat"
)

colnames(coords_test) <- c(
    "Lon",
    "Lat"
)


# ============================================================
# 5. Scale spatial coordinates to [0,1]^2
#
# IMPORTANT:
# Training range is used for BOTH training and testing.
# ============================================================

lon_min <- min(coords_train[, 1])
lon_max <- max(coords_train[, 1])

lat_min <- min(coords_train[, 2])
lat_max <- max(coords_train[, 2])


coords_train_scaled <- coords_train

coords_train_scaled[, 1] <-
    (coords_train[, 1] - lon_min) /
    (lon_max - lon_min)

coords_train_scaled[, 2] <-
    (coords_train[, 2] - lat_min) /
    (lat_max - lat_min)


coords_test_scaled <- coords_test

coords_test_scaled[, 1] <-
    (coords_test[, 1] - lon_min) /
    (lon_max - lon_min)

coords_test_scaled[, 2] <-
    (coords_test[, 2] - lat_min) /
    (lat_max - lat_min)


# ============================================================
# 6. Standardize temperature
#
# IMPORTANT:
# mean and SD are calculated ONLY from training data.
#
# y_standardized = (y - mean_train) / sd_train
# ============================================================

temp_mean <- mean(
    D2$MaskTemp
)

temp_sd <- sd(
    D2$MaskTemp
)


y_train <- (
    D2$MaskTemp - temp_mean
) / temp_sd


y_train <- matrix(
    y_train,
    ncol = 1
)


cat("\nTemperature standardization:\n")
cat("Training mean =", temp_mean, "\n")
cat("Training SD   =", temp_sd, "\n")


# ============================================================
# 7. 2D Wendland basis function
# ============================================================

wendland_basis <- function(
    coords,
    resolutions = c(10, 19, 37, 73),
    support_factor = 2.5
) {

    N <- nrow(coords)

    basis_list <- list()


    for (res in resolutions) {

        # ----------------------------------------------------
        # 2D regular grid of knots
        # ----------------------------------------------------

        x_knots <- seq(
            0,
            1,
            length.out = res
        )

        y_knots <- seq(
            0,
            1,
            length.out = res
        )


        knots <- expand.grid(
            x = x_knots,
            y = y_knots
        )

        knots <- as.matrix(knots)


        # ----------------------------------------------------
        # Knot spacing
        # ----------------------------------------------------

        spacing <- 1 / (res - 1)


        # ----------------------------------------------------
        # Support radius
        # ----------------------------------------------------

        theta <- support_factor * spacing


        # ----------------------------------------------------
        # Euclidean distance
        # ----------------------------------------------------

        dx <- outer(
            coords[, 1],
            knots[, 1],
            "-"
        )

        dy <- outer(
            coords[, 2],
            knots[, 2],
            "-"
        )


        distance <- sqrt(
            dx^2 + dy^2
        )


        # ----------------------------------------------------
        # Scaled distance
        # ----------------------------------------------------

        r <- distance / theta


        # ----------------------------------------------------
        # Wendland function
        # ----------------------------------------------------

        phi <- matrix(
            0,
            nrow = N,
            ncol = nrow(knots)
        )


        ind <- r <= 1


        phi[ind] <-
            (1 - r[ind])^6 *
            (
                35 * r[ind]^2 +
                18 * r[ind] +
                3
            ) / 3


        basis_list[[length(basis_list) + 1]] <-
            phi


        cat(
            "Resolution:",
            res, "x", res,
            "| Basis:",
            res^2,
            "| theta:",
            theta,
            "\n"
        )
    }


    # --------------------------------------------------------
    # Combine all resolutions
    # --------------------------------------------------------

    phi_all <- do.call(
        cbind,
        basis_list
    )


    return(phi_all)
}


# ============================================================
# 8. Create training basis
# ============================================================

resolutions <- c(
    10,
    19,
    37,
    73
)


phi_train <- wendland_basis(
    coords = coords_train_scaled,
    resolutions = resolutions,
    support_factor = 2.5
)


# ============================================================
# 9. Create testing basis
#
# Same resolutions and same basis construction.
# ============================================================

phi_test <- wendland_basis(
    coords = coords_test_scaled,
    resolutions = resolutions,
    support_factor = 2.5
)

t2 <- proc.time()-t1
t2

# ============================================================
# 10. Check dimensions
# ============================================================

K <- ncol(phi_train)


cat("\n============================================\n")
cat("Training observations :", nrow(phi_train), "\n")
cat("Training basis        :", ncol(phi_train), "\n")
cat("Testing observations  :", nrow(phi_test), "\n")
cat("Testing basis         :", ncol(phi_test), "\n")
cat("============================================\n\n")


stopifnot(
    ncol(phi_train) == ncol(phi_test)
)


# ============================================================
# 11. DeepKriging neural network
#
# Three hidden layers
# Each hidden layer = 100 neurons
# ============================================================

t1 <- proc.time()

model <- keras_model_sequential() %>%

    # Hidden layer 1
    layer_dense(
        units = 100,
        activation = "relu",
        input_shape = c(K)
    ) %>%

    # Hidden layer 2
    layer_dense(
        units = 100,
        activation = "relu"
    ) %>%

    # Hidden layer 3
    layer_dense(
        units = 100,
        activation = "relu"
    ) %>%

    # Output layer
    layer_dense(
        units = 1,
        activation = "linear"
    )


# ============================================================
# 12. Adam optimizer
#
# Learning rate = 0.001
# ============================================================

optimizer <- optimizer_adam(
    learning_rate = 0.001
)


model %>%
    compile(
        loss = "mse",
        optimizer = optimizer,
        metrics = c("mse")
    )


# ============================================================
# 13. Early stopping
#
# Stop if validation loss does not improve for 20 epochs.
# ============================================================

early_stop <- callback_early_stopping(
    monitor = "val_loss",
    patience = 20,
    restore_best_weights = TRUE
)


# ============================================================
# 14. Train
# ============================================================

history <- model %>%

    fit(
        phi_train,
        y_train,

        epochs = 200,

        batch_size = 1024,

        validation_split = 0.1,

        callbacks = list(
            early_stop
        ),

        verbose = 1
    )

t3 <- proc.time()-t1
t3

# ============================================================
# 15. Training prediction
#
# Prediction is on standardized scale.
# Then transform back to original temperature scale.
# ============================================================

pred_train_std <- model %>%
    predict(phi_train)


pred_train <- 
    as.vector(pred_train_std) *
    temp_sd +
    temp_mean


# ============================================================
# 16. Testing prediction
# ============================================================

t1 <- proc.time()

pred_test_std <- model %>%
    predict(phi_test)


pred_test <-
    as.vector(pred_test_std) *
    temp_sd +
    temp_mean

t4 <- proc.time()-t1
t4

# ============================================================
# 17. Test truth
# ============================================================

y_test <- D3$TrueTemp


# ============================================================
# 18. Calculate MAE and RMSE
# ============================================================

# ------------------------------------------------------------
# Training
# ------------------------------------------------------------

MAE_train <- mean(
    abs(
        D2$MaskTemp -
        pred_train
    )
)


RMSE_train <- sqrt(
    mean(
        (
            D2$MaskTemp -
            pred_train
        )^2
    )
)


# ------------------------------------------------------------
# Testing
# ------------------------------------------------------------

MAE_test <- mean(
    abs(
        y_test -
        pred_test
    )
)


RMSE_test <- sqrt(
    mean(
        (
            y_test -
            pred_test
        )^2
    )
)


# ============================================================
# 19. Print prediction performance
# ============================================================

cat("\n============================================\n")
cat("DeepKriging Prediction Performance\n")
cat("============================================\n")

cat("\nTraining:\n")
cat("MAE  =", MAE_train, "\n")
cat("RMSE =", RMSE_train, "\n")

cat("\nTesting:\n")
cat("MAE  =", MAE_test, "\n")
cat("RMSE =", RMSE_test, "\n")

cat("\n============================================\n")


# ============================================================
# 20. Create output directory
# ============================================================

output_dir <- "DeepKriging_Python"


if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}


# ============================================================
# 21. Export files for Python
# ============================================================


# ------------------------------------------------------------
# Training response
# ------------------------------------------------------------

write.csv(
    data.frame(
        y_train = D2$MaskTemp
    ),
    file.path(
        output_dir,
        "y_train.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------
# Training coordinates
# ------------------------------------------------------------

write.csv(
    data.frame(
        Lon = coords_train[, 1],
        Lat = coords_train[, 2]
    ),
    file.path(
        output_dir,
        "coords_train.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------
# Training prediction
# ------------------------------------------------------------

write.csv(
    data.frame(
        pred_train = pred_train
    ),
    file.path(
        output_dir,
        "pred_train.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------
# Testing coordinates
# ------------------------------------------------------------

write.csv(
    data.frame(
        Lon = coords_test[, 1],
        Lat = coords_test[, 2]
    ),
    file.path(
        output_dir,
        "coords_test.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------
# Testing prediction
# ------------------------------------------------------------

write.csv(
    data.frame(
        pred_test = pred_test
    ),
    file.path(
        output_dir,
        "pred_test.csv"
    ),
    row.names = FALSE
)


# ------------------------------------------------------------
# Testing TRUE values
# ------------------------------------------------------------

write.csv(
    data.frame(
        y_test = y_test
    ),
    file.path(
        output_dir,
        "y_test.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# Wendland basis for DCDR
# ------------------------------------------------------------

write.csv(
    phi_train,
    file.path(
        output_dir,
        "phi_train.csv"
    ),
    row.names = FALSE
)

write.csv(
    phi_test,
    file.path(
        output_dir,
        "phi_test.csv"
    ),
    row.names = FALSE
)


# ============================================================
# 22. Print output information
# ============================================================

cat("\n============================================\n")
cat("Files exported to:", output_dir, "\n")
cat("============================================\n")

cat(
    "y_train.csv      :",
    nrow(D2),
    "rows\n"
)

cat(
    "coords_train.csv :",
    nrow(coords_train),
    "rows\n"
)

cat(
    "pred_train.csv   :",
    length(pred_train),
    "rows\n"
)

cat(
    "coords_test.csv  :",
    nrow(coords_test),
    "rows\n"
)

cat(
    "pred_test.csv    :",
    length(pred_test),
    "rows\n"
)

cat(
    "y_test.csv       :",
    length(y_test),
    "rows\n"
)

cat("\nCPU cores =", num_cores, "\n")

cat("\n============================================\n")


# ============================================================
# 23. Plot training history
# ============================================================

plot(history)

save.image("DeepKriging.RData")