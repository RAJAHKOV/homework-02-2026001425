Rcpp::sourceCpp("/home/rstudio/homework2/Kalman-smoothing/Kalman-filtering.cpp")
Rcpp::sourceCpp("/home/rstudio/homework2/Kalman-smoothing/Kalman-smoothing.cpp")
start = Sys.time()
generateData <- function(n, A, C, z0, sigmaQ, sigmaR) {
    z <- z0
    ZActual <- matrix(0, length(z0), n)
    Yobs <- matrix(0, nrow(C), n)
    for (i in seq_len(n)) {
        z <- as.vector(A %*% z) + rnorm(length(z0), sd = sigmaQ)
        Yobs[, i] <- as.vector(C %*% z) + rnorm(nrow(C), sd = sigmaR)
        ZActual[, i] <- z
    }
    list(Yobs = Yobs, ZActual = ZActual)
}

set.seed(2026)
sigmaQ <- 0.2
sigmaR <- 0.8
n <- 100L
A <- diag(4)
A[1, 3] <- A[2, 4] <- 0.5
C <- diag(1, nrow = 2, ncol = 4)
z0 <- c(0, 0, 1, 1)
kdata <- generateData(n, A, C, z0, sigmaQ, sigmaR)
filtered <- kalmanFiltering(kdata$Yobs, z0, A, C, sigmaQ, sigmaR)
smoothed <- kalmanSmoothing(kdata$Yobs, z0, A, C, sigmaQ, sigmaR)
stopifnot(identical(dim(smoothed), c(4L, n)), all(is.finite(smoothed)),
          max(abs(smoothed[, n] - filtered[, n])) < 1e-10)

# One observation: filtering and smoothing should coincide.
single <- kalmanSmoothing(kdata$Yobs[, 1, drop = FALSE], z0, A, C, sigmaQ, sigmaR)
stopifnot(max(abs(single[, 1] - filtered[, 1])) < 1e-10)
print(c(filtering_mse = mean((filtered[1:2, ] - kdata$ZActual[1:2, ])^2),
        smoothing_mse = mean((smoothed[1:2, ] - kdata$ZActual[1:2, ])^2)))

end <- Sys.time()
cat("kalman-smoothing 计算耗时",end-start,'s')


drawPath <- function(estimate, title) {
    plot(t(kdata$Yobs), pch = 20, col = "grey70", xlab = "x", ylab = "y",
         main = title, 
         #asp = 1,
         #xlim = range(kdata$Yobs[1, ], kdata$ZActual[1, ], estimate[1, ]),
         #ylim = range(kdata$Yobs[2, ], kdata$ZActual[2, ], estimate[2, ])
         xlim = c(-20,20),   
         ylim = c(-5,120)
    )
    lines(t(kdata$ZActual[1:2, ]), col = "black", lwd = 2)
    lines(t(estimate[1:2, ]), col = "blue", lwd = 2)
    legend("topleft", c("Observed", "True", "Estimated"),
           col = c("grey70", "black", "blue"), pch = c(20, NA, NA),
           lty = c(NA, 1, 1), bty = "n")
}
old_par <- par(mfrow = c(2, 1), mar = c(3, 3, 3, 1))
drawPath(filtered, "Kalman filtering")
drawPath(smoothed, "Kalman smoothing")
par(old_par)
