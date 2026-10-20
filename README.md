# homework-02-2026001425

## 1.image-convolution
```bash
source("~/homework2/convolution/image-convolution.R")
```
返回
```bash
$gaussian_differences
[1] 0

$sharpen_differences
[1] 0

$edge_differences
[1] 0

  method seconds repetitions rows cols
1    C++  0.0023          10  256  256
2      R  0.4888          10  256  256
```
## 2.log-mvn-density
```bash
source("~/homework2/mvn_log_density/mvn-log-density.R")
```
返回
```bash
actual  expected 
-1.504404 -1.504404 
mvn-log-likelihood计算耗时： 0.0005311966 s
   actual  expected 
-5.079362 -5.079362 
mvn-log-likelihood计算耗时： 0.0009565353 s
   actual  expected 
-4.096459 -4.096459 
mvn-log-likelihood计算耗时： 0.0001516342 s
   actual  expected 
-230.0917 -230.0917 
mvn-log-likelihood计算耗时： 0.0001401901 s
```

## 3.Kalman Smoothing
```bash
source("~/homework2/Kalman-smoothing/kalman-smoothing.R")
```
返回
```bash
filtering_mse smoothing_mse 
   0.25161834    0.06381335 
kalman-smoothing 计算耗时 0.007317543 s
```

## 4.Ridge Regression and Tikhonov Regularization
```bash
source("~/homework2/Ridge Regression and Tikhonov Regression/ridge-cv.R")
```
返回
```bash
direct 
actual =  4.196365 
expected =  4.196365 
diff =  0 
leverage 
actual =  4.196365 
expected =  4.196365 
diff =  2.664535e-15 
svd 
actual =  4.196365 
expected =  4.196365 
diff =  8.881784e-16 
           method max_absolute_error selected_lambda
direct     direct       0.000000e+00       0.6951928
leverage leverage       1.421085e-13       0.6951928
svd           svd       1.705303e-12       0.6951928
           method median_seconds repetitions calls_per_batch   n  p lambda_count
direct     direct       0.060000           5               1 200 20           20
leverage leverage       0.231170           5             100 200 20           20
svd           svd       0.000176           5             500 200 20           20
```

```bash
source("~/homework2/Ridge Regression and Tikhonov Regularization/tikhonov.R")
```
返回
```bash
coefficient_error      fitted_error consistency_error 
     2.220446e-15      7.105427e-15      6.217249e-15 
ridge_max_absolute_error 
            4.440892e-15
```

5.AI使用情况
在四道题中均有使用ChatGPT-6，主要用于改正自己写的代码的错误，以及文件管理错误。通过比较运行结果的差异和无Warning提醒，我验证结果是正确的。
