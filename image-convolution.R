Rcpp::sourceCpp("~/homework2/convolution/image-convolution.cpp")

# R uses indices 1, ..., length, unlike the C++ helper.
mirrorIndexR <- function(index, length) {
  if (index < 1L) return(1L - index)
  if (index > length) return(2L * length + 1L - index)
  index
}

imageConvR <- function(img, kernel) {
  # TODO: Use loops with the same boundary and clipping rules as imageConv.
  H <- nrow(img)
  W <- ncol(img)
  out <- matrix(0.0, nrow = H, ncol = W)
  
  # compute every element of out[].note that the coordinates of Vectors in R starts from 1
  for (i in 1:H) {
    for (j in 1:W) {
      sum_val <- 0.0
      for (a in -1:1) {
        for (b in -1:1) {
          # image coordinates 
          r <- mirrorIndexR(i + a, H)
          c <- mirrorIndexR(j + b, W)
          # The subscripts of the kernel are a+2, b+2,ranging 1,2,3
          k_val <- kernel[a + 2, b + 2]
          sum_val <- sum_val + img[r, c] * k_val
        }
      }
      # sum_val is truncated into [0,1]
      if (sum_val < 0) sum_val <- 0
      if (sum_val > 1) sum_val <- 1
      out[i,j] <- sum_val
    }
  }
  return(out)
}

# A reproducible grayscale test image, with no extra image-package dependency.
grid <- seq(-1, 1, length.out = 256)
img <- outer(grid, grid, function(x, y) {
  0.2 + 0.5 * (x^2 + y^2 < 0.5) + 0.2 * (abs(x) < 0.2)
})
kernels <- list(
  gaussian = matrix(c(1, 2, 1, 2, 4, 2, 1, 2, 1), 3, 3) / 16,
  sharpen = matrix(c(0, -1, 0, -1, 5, -1, 0, -1, 0), 3, 3),
  edge = matrix(c(-1, -1, -1, -1, 8, -1, -1, -1, -1), 3, 3)
)


results <- lapply(kernels, function(kernel) {
  cpp <- imageConv(img, kernel)
  r <- imageConvR(img, kernel)
  stopifnot(identical(dim(cpp), dim(img)),
            all(cpp >= 0 & cpp <= 1), max(abs(cpp - r)) < 1e-10)
  cpp
})

ans <- vapply(seq_along(kernels), function(i) {
  max(abs(results[[i]] - imageConvR(img, kernels[[i]])))
}, numeric(1))

ans <- list(
  gaussian_differences = ans[1],
  sharpen_differences = ans[2],
  edge_differences = ans[3]
)

print(ans)

#cat("cpp",dim(results[[1]]),'\n')
#cat("R:",dim(imageConvR(img,kernels[[1]])))


# Compile and warm up outside timing; report elapsed seconds per call.
repetitions <- 10L
kernel <- kernels$gaussian
invisible(imageConv(img, kernel))
invisible(imageConvR(img, kernel))
cpp_time <- system.time(for (i in seq_len(repetitions)) {
  imageConv(img, kernel)
})[["elapsed"]] / repetitions
r_time <- system.time(for (i in seq_len(repetitions)) {
  imageConvR(img, kernel)
})[["elapsed"]] / repetitions
print(data.frame(method = c("C++", "R"), seconds = c(cpp_time, r_time),
                 repetitions = repetitions, rows = nrow(img), cols = ncol(img)))

drawImage <- function(values, title) {
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), asp = 1)
  rasterImage(as.raster(values), 0, 0, 1, 1, interpolate = FALSE)
  title(main = title)
}
old_par <- par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))
drawImage(img, "Original")
for (name in names(results)) drawImage(results[[name]], name)
par(old_par)