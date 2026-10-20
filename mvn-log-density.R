Rcpp::sourceCpp("/home/rstudio/homework2/mvn_log_density/mvn-log-density.cpp")

checkDensity <- function(x, mu, Sigma, expected) {
  start <- Sys.time()
  actual <- mvnLogDensity(x, mu, Sigma)
  stopifnot(length(actual) == 1L, is.finite(actual),
            abs(actual - expected) < 1e-10)
  print(c(actual = actual, expected = expected))
  end <- Sys.time()
  cat("mvn-log-likelihood计算耗时：", end - start, "s\n")
}

# One dimension: compare with base R's normal log-density.
checkDensity(1.2, 0.3, matrix(2.25, 1, 1),
             dnorm(1.2, mean = 0.3, sd = 1.5, log = TRUE))

# Independent coordinates: sum univariate log-densities.
x <- c(1, -2, 0.5)
mu <- c(0.2, 0, -0.3)
variances <- c(1, 2, 3)
checkDensity(x, mu, diag(variances),
             sum(dnorm(x, mu, sqrt(variances), log = TRUE)))

# A non-diagonal covariance: construct Sigma = L L^T with known L.
L <- matrix(c(1.5, 0.4, -0.2, 0, 1.2, 0.3, 0, 0, 0.8), 3, 3)
v <- c(0.5, -0.7, 1.1)
Sigma <- L %*% t(L)
x <- as.vector(mu + L %*% v)
expected <- -length(x) / 2 * log(2 * pi) - sum(log(diag(L))) - sum(v^2) / 2
checkDensity(x, mu, Sigma, expected)

# TODO: Add at least one test of your own.
# An ordinary case
S <- matrix(c(6.0,4.3,4.3,9.0),nrow=2,ncol=2)
M <- c(44.4,67.8)
x <- c(0.4,8.5)
expected <- -1*log(2*pi)-0.5*log(6.0*9.0-4.3*4.3)-0.5*t(x-M)%*%solve(S)%*%(x-M)
checkDensity(x, M, S, expected)
