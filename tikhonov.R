Rcpp::sourceCpp("~/homework2/Ridge Regression and Tikhonov Regularization/tikhonov.cpp")

# computeGSVD() is already implemented in tikhonov.cpp and exported to R.
# Only tikhonovPath() needs to be completed by the student.

checkFactors <- function(X, B, factors) {
    p <- ncol(X)
    stopifnot(max(abs(crossprod(factors$U) - diag(p))) < 1e-10,
              max(abs(crossprod(factors$U2) - diag(p))) < 1e-10,
              max(abs(X %*% factors$Z - factors$U %*%
                      diag(factors$alpha, p, p))) < 1e-10,
              max(abs(B %*% factors$Z - factors$U2 %*%
                      diag(factors$beta, p, p))) < 1e-10,
              max(abs(factors$alpha^2 + factors$beta^2 - 1)) < 1e-10)
}

set.seed(2026)
n <- 100L
p <- 6L
X <- matrix(rnorm(n * p), n, p)
B <- diag(seq(1, 2, length.out = p))
B[cbind(seq_len(p - 1L), 2:p)] <- -0.3
beta_true <- c(1, -2, 0.5, 0, 1.5, -0.5)
y <- as.vector(X %*% beta_true + rnorm(n, sd = 0.5))
lambdas <- c(0.01, 0.1, 1, 10, 100)

factors <- computeGSVD(X, B)
checkFactors(X, B, factors)
result <- tikhonovPath(factors$U, factors$Z, factors$alpha, factors$beta,
                       y, lambdas)
reference <- vapply(lambdas, function(lambda) {
    as.vector(solve(crossprod(X) + lambda * crossprod(B), crossprod(X, y)))
}, numeric(p))
stopifnot(is.list(result), all(c("coefficients", "fitted") %in% names(result)),
          identical(dim(result$coefficients), c(p, length(lambdas))),
          identical(dim(result$fitted), c(n, length(lambdas))),
          all(is.finite(result$coefficients)), all(is.finite(result$fitted)))
errors <- c(coefficient_error = max(abs(result$coefficients - reference)),
            fitted_error = max(abs(result$fitted - X %*% reference)),
            consistency_error = max(abs(result$fitted - X %*% result$coefficients)))
print(errors)
stopifnot(max(errors) < 1e-8)

# Compute factors for B = I with the same supplied helper.
ridge_factors <- computeGSVD(X, diag(p))
checkFactors(X, diag(p), ridge_factors)
ridge <- tikhonovPath(ridge_factors$U, ridge_factors$Z,
                      ridge_factors$alpha, ridge_factors$beta, y, lambdas)
# Independently compare with ordinary ridge regression via an SVD of X.
decomposition <- svd(X)
ridge_reference <- vapply(lambdas, function(lambda) {
    as.vector(decomposition$v %*%
              ((decomposition$d / (decomposition$d^2 + lambda)) *
               as.vector(crossprod(decomposition$u, y))))
}, numeric(p))
stopifnot(identical(dim(ridge$coefficients), c(p, length(lambdas))),
          identical(dim(ridge$fitted), c(n, length(lambdas))),
          all(is.finite(ridge$coefficients)), all(is.finite(ridge$fitted)))
ridge_error <- max(abs(ridge$coefficients - ridge_reference),
                   abs(ridge$fitted - X %*% ridge_reference))
print(c(ridge_max_absolute_error = ridge_error))
stopifnot(ridge_error < 1e-8)
