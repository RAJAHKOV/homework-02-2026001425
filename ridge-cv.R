Rcpp::sourceCpp("~/homework2/Ridge Regression and Tikhonov Regularization/ridge-cv.cpp")

methods <- list(direct = ridgeCVDirect, leverage = ridgeCVLeverage, svd = ridgeCVSVD)

# A small test with an independently calculated LOOCV value.
small_X <- matrix(c(1, 2, 3), ncol = 1)
small_y <- c(1, 0, 2)
small_expected <- (4 / 7)^2 + (14 / 11)^2 + (3 / 2)^2

for (method in methods) {
    actual <- as.vector(method(small_X, small_y, 1))
    
    stopifnot(length(actual) == 1L, is.finite(actual),
             abs(actual - small_expected) < 1e-10)
    
}

set.seed(2026)
n <- 200L
p <- 20L
X <- matrix(rnorm(n * p), n, p)
beta_true <- c(1, -2, 0.5, 1.5, -0.5, rep(0, p - 5L))
y <- as.vector(X %*% beta_true + rnorm(n, sd = 0.5))
lambdas <- 10^seq(-3, 3, length.out = 20)

cv <- vapply(methods, function(method) {
    values <- as.vector(method(X, y, lambdas))
    stopifnot(length(values) == length(lambdas), all(is.finite(values)),
              all(values >= 0))
    values
}, numeric(length(lambdas)))
differences <- apply(cv, 2, function(values) max(abs(values - cv[, "direct"])))
print(data.frame(method = names(methods), max_absolute_error = differences,
                 selected_lambda = lambdas[apply(cv, 2, which.min)]))
stopifnot(max(differences) < 1e-8)

# Each timed call includes all required preprocessing and factorization.
# Compile, generate data, and warm up before timing.
for (method in methods) invisible(method(X, y, lambdas))
repetitions <- 5L
# Batch fast methods to avoid rounding very short runs to zero seconds.
# Divide each batch time by its call count for a per-call comparison.
calls_per_batch <- c(direct = 1L, leverage = 100L, svd = 500L)
timing <- matrix(NA_real_, repetitions, length(methods),
                 dimnames = list(NULL, names(methods)))
for (i in seq_len(repetitions)) {
    for (j in sample(seq_along(methods))) {
        batch_calls <- calls_per_batch[j]
        timing[i, j] <- system.time(for (k in seq_len(batch_calls)) {
            invisible(methods[[j]](X, y, lambdas))
        })[["elapsed"]] / batch_calls
    }
}
print(data.frame(method = names(methods),
                 median_seconds = apply(timing, 2, median),
                 repetitions = repetitions, calls_per_batch = calls_per_batch,
                 n = n, p = p, lambda_count = length(lambdas)))

matplot(lambdas, cv, type = "b", log = "x", lty = 1:3, pch = 1:3,
        col = c("black", "blue", "red"), xlab = "lambda", ylab = "LOOCV error")
legend("topleft", names(methods), col = c("black", "blue", "red"),
       lty = 1:3, pch = 1:3, bty = "n")
