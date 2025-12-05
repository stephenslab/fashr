
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fashr

<!-- badges: start -->

<!-- badges: end -->

The goal of fashr is to implement the functional adaptive shrinkage
through empirical Bayes.

**Note:** The code to generate the results for the manuscript is
[here](https://github.com/stephenslab/fashr-paper).

## Installation

You can install the development version of fashr from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("stephenslab/fashr")
#> 
#> ✔ Updated metadata database: 7.63 MB in 10 files.
#> 
#> ℹ Updating metadata database
#> ✔ Updating metadata database ... done
#> 
#> 
#> → Will install 28 packages.
#> → Will update 1 package.
#> → Will download 27 CRAN packages (45.39 MB), cached: 1 (4.10 MB).
#> → Will download 1 package with unknown size.
#> + cli                    3.6.5      🔧 ⬇ (1.47 MB)
#> + cowplot                1.2.0       ⬇ (1.38 MB)
#> + farver                 2.1.2      🔧 ⬇ (1.97 MB)
#> + fashr         0.1.39 → 0.1.39     👷🏻‍♂️🔧 ⬇ (GitHub: 3c3390b)
#> + ggplot2                4.0.1       ⬇ (8.47 MB)
#> + glue                   1.8.0      🔧 ⬇ (175.14 kB)
#> + gtable                 0.3.6       ⬇ (225.12 kB)
#> + irlba                  2.3.5.1    🔧 ⬇ (319.80 kB)
#> + isoband                0.2.7      🔧 ⬇ (1.87 MB)
#> + labeling               0.4.3       ⬇ (61.83 kB)
#> + LaplacesDemon          16.1.6     
#> + lifecycle              1.0.4       ⬇ (125.64 kB)
#> + magrittr               2.0.4      🔧 ⬇ (231.97 kB)
#> + mixsqp                 0.3-54     🔧 ⬇ (994.47 kB)
#> + numDeriv               2016.8-1.1  ⬇ (115.06 kB)
#> + plyr                   1.8.9      🔧 ⬇ (984.75 kB)
#> + R6                     2.6.1       ⬇ (87.28 kB)
#> + RColorBrewer           1.1-3       ⬇ (51.79 kB)
#> + Rcpp                   1.1.0      🔧 ⬇ (3.37 MB)
#> + reshape2               1.4.5      🔧 ⬇ (343.07 kB)
#> + rlang                  1.1.6      🔧 ⬇ (1.91 MB)
#> + S7                     0.2.1      🔧 ⬇ (343.27 kB)
#> + scales                 1.4.0       ⬇ (873.34 kB)
#> + stringi                1.8.7      🔧 ⬇ (14.79 MB)
#> + stringr                1.6.0       ⬇ (333.11 kB)
#> + TMB                    1.9.18     🔧 ⬇ (1.47 MB)
#> + vctrs                  0.6.5      🔧 ⬇ (1.89 MB)
#> + viridisLite            0.4.2       ⬇ (1.30 MB)
#> + withr                  3.0.2       ⬇ (224.91 kB)
#> ℹ Getting 27 pkgs (45.39 MB) and 1 pkg with unknown size, 1 (4.10 MB) cached
#> ✔ Got R6 2.6.1 (aarch64-apple-darwin20) (87.28 kB)
#> ✔ Got S7 0.2.1 (aarch64-apple-darwin20) (343.27 kB)
#> ✔ Got RColorBrewer 1.1-3 (aarch64-apple-darwin20) (51.79 kB)
#> ✔ Got labeling 0.4.3 (aarch64-apple-darwin20) (61.83 kB)
#> ✔ Got gtable 0.3.6 (aarch64-apple-darwin20) (225.12 kB)
#> ✔ Got cli 3.6.5 (aarch64-apple-darwin20) (1.47 MB)
#> ✔ Got magrittr 2.0.4 (aarch64-apple-darwin20) (231.97 kB)
#> ✔ Got withr 3.0.2 (aarch64-apple-darwin20) (224.91 kB)
#> ✔ Got TMB 1.9.18 (aarch64-apple-darwin20) (1.47 MB)
#> ✔ Got Rcpp 1.1.0 (aarch64-apple-darwin20) (3.37 MB)
#> ✔ Got vctrs 0.6.5 (aarch64-apple-darwin20) (1.89 MB)
#> ✔ Got lifecycle 1.0.4 (aarch64-apple-darwin20) (125.64 kB)
#> ✔ Got numDeriv 2016.8-1.1 (aarch64-apple-darwin20) (115.06 kB)
#> ✔ Got stringr 1.6.0 (aarch64-apple-darwin20) (333.11 kB)
#> ✔ Got cowplot 1.2.0 (aarch64-apple-darwin20) (1.38 MB)
#> ✔ Got plyr 1.8.9 (aarch64-apple-darwin20) (984.75 kB)
#> ✔ Got irlba 2.3.5.1 (aarch64-apple-darwin20) (319.80 kB)
#> ✔ Got glue 1.8.0 (aarch64-apple-darwin20) (175.14 kB)
#> ✔ Got viridisLite 0.4.2 (aarch64-apple-darwin20) (1.30 MB)
#> ✔ Got isoband 0.2.7 (aarch64-apple-darwin20) (1.87 MB)
#> ✔ Got mixsqp 0.3-54 (aarch64-apple-darwin20) (994.47 kB)
#> ✔ Got reshape2 1.4.5 (aarch64-apple-darwin20) (343.07 kB)
#> ✔ Got rlang 1.1.6 (aarch64-apple-darwin20) (1.91 MB)
#> ✔ Got farver 2.1.2 (aarch64-apple-darwin20) (1.97 MB)
#> ✔ Got scales 1.4.0 (aarch64-apple-darwin20) (873.34 kB)
#> ✔ Got stringi 1.8.7 (aarch64-apple-darwin20) (14.79 MB)
#> ✔ Got fashr 0.1.39 (source) (3.07 MB)
#> ✔ Got ggplot2 4.0.1 (aarch64-apple-darwin20) (8.47 MB)
#> ✔ Installed LaplacesDemon 16.1.6  (143ms)
#> ✔ Installed R6 2.6.1  (147ms)
#> ✔ Installed RColorBrewer 1.1-3  (149ms)
#> ✔ Installed S7 0.2.1  (139ms)
#> ✔ Installed cli 3.6.5  (121ms)
#> ✔ Installed cowplot 1.2.0  (120ms)
#> ✔ Installed Rcpp 1.1.0  (213ms)
#> ✔ Installed TMB 1.9.18  (224ms)
#> ✔ Installed farver 2.1.2  (190ms)
#> ✔ Installed ggplot2 4.0.1  (193ms)
#> ✔ Installed glue 1.8.0  (112ms)
#> ✔ Installed gtable 0.3.6  (42ms)
#> ✔ Installed irlba 2.3.5.1  (41ms)
#> ✔ Installed isoband 0.2.7  (42ms)
#> ✔ Installed labeling 0.4.3  (40ms)
#> ✔ Installed lifecycle 1.0.4  (40ms)
#> ✔ Installed magrittr 2.0.4  (41ms)
#> ✔ Installed mixsqp 0.3-54  (44ms)
#> ✔ Installed numDeriv 2016.8-1.1  (67ms)
#> ✔ Installed plyr 1.8.9  (43ms)
#> ✔ Installed reshape2 1.4.5  (43ms)
#> ✔ Installed rlang 1.1.6  (43ms)
#> ✔ Installed scales 1.4.0  (44ms)
#> ✔ Installed stringr 1.6.0  (22ms)
#> ✔ Installed stringi 1.8.7  (92ms)
#> ✔ Installed vctrs 0.6.5  (43ms)
#> ✔ Installed viridisLite 0.4.2  (42ms)
#> ✔ Installed withr 3.0.2  (25ms)
#> ℹ Packaging fashr 0.1.39
#> ✔ Packaged fashr 0.1.39 (714ms)
#> ℹ Building fashr 0.1.39
#> ✔ Built fashr 0.1.39 (1m 37.9s)
#> ✔ Installed fashr 0.1.39 (github::stephenslab/fashr@3c3390b) (80ms)
#> ✔ 1 pkg + 30 deps: kept 1, upd 1, added 28, dld 28 (NA B) [1m 51.7s]
```

## Example

This is a basic example, which tries to identify non-linear functions in
a set of datasets. We use `order = 2` to apply the second order IWP,
which has a base model S<sub>0</sub> being the space of linear
functions.

``` r
library(fashr)

# Simulate five datasets
set.seed(123)
datasets <- list()
for (i in 1:20) {
  n <- 50
  t <- seq(0, 5, length.out = n)
  sd <- sample(c(2, 1), size = n, replace = TRUE)
  u <- runif(1); if (u < 0.5) {f <- function(t) 3*cos(t)} else {f <- function(t) (t)}
  y <- f(t) + rnorm(n, sd = 0.5)
  datasets[[i]] <- data.frame(t = t, y = y, sd = sd)
}

# Fit the model
fash_fit <- fash(Y = "y", smooth_var = "t", S = "sd", data = datasets, 
                  order = 2, likelihood = "gaussian", verbose = TRUE)
```

``` r
fash_fit
#> Fitted fash Object
#> -------------------
#> Number of datasets: 20
#> Likelihood: gaussian
#> Number of PSD grid values: 25 (initial), 2 (non-trivial)
#> Order of Integrated Wiener Process (IWP): 2
```

Take a look at the structure plot ordered by local false discovery rate:

``` r
plot(fash_fit, ordering = "lfdr")
```

<img src="man/figures/README-structure_plot-1.png" width="100%" />

Obtain the posterior summary of the function for the first dataset:

``` r
fitted <- predict(fash_fit, index = 1)
str(fitted)
#> 'data.frame':    50 obs. of  5 variables:
#>  $ x     : num  0 0.102 0.204 0.306 0.408 ...
#>  $ mean  : num  3.43 3.2 2.98 2.75 2.52 ...
#>  $ median: num  3.43 3.2 2.98 2.75 2.53 ...
#>  $ lower : num  1.96 1.94 1.89 1.8 1.67 ...
#>  $ upper : num  4.88 4.44 4.04 3.68 3.36 ...
```

Obtain the posterior samples of the function:

``` r
fitted_samps <- predict(fash_fit, index = 1, only.samples = TRUE, M = 30)
str(fitted_samps)
#>  num [1:50, 1:30] 3.36 3.13 2.96 2.78 2.57 ...
```

Visualize the predictions for the first dataset together with the data:

``` r
plot(datasets[[1]]$t, datasets[[1]]$y, type = "p", col = "black", ylab = "y", xlab = "t")
lines(fitted$x, fitted$mean, col = "red")
polygon(c(fitted$x, rev(fitted$x)), c(fitted$lower, rev(fitted$upper)), col = rgb(1, 0, 0, 0.2), border = NA)
matlines(fitted$x, fitted_samps[, 1:5], col = "blue", lty = 2, lwd = 0.5)
```

<img src="man/figures/README-predictions-1.png" width="100%" />

Compute the FDR, and highlight the datasets with FDR \< 0.1:

``` r
fdr_result <- fdr_control(fash_fit, alpha = 0.1, plot = TRUE)
#> 5 datasets are significant at alpha level 0.10. Total datasets tested: 20.
```

<img src="man/figures/README-fdr-1.png" width="100%" />

``` r
fdr_result
#> $fdr_results
#>            rank index         lfdr          FDR
#> Dataset_1     5     1 1.981230e-08 3.988248e-09
#> Dataset_2    18     2 9.794739e-01 7.023167e-01
#> Dataset_3     8     3 9.683893e-01 3.609240e-01
#> Dataset_4     1     4 3.492031e-16 3.492031e-16
#> Dataset_5    11     5 9.719933e-01 5.275334e-01
#> Dataset_6    15     6 9.778700e-01 6.470013e-01
#> Dataset_7    19     7 9.803864e-01 7.169520e-01
#> Dataset_8    17     8 9.788706e-01 6.860134e-01
#> Dataset_9     7     9 9.673522e-01 2.741432e-01
#> Dataset_10   14    10 9.752621e-01 6.233678e-01
#> Dataset_11   13    11 9.748583e-01 5.962990e-01
#> Dataset_12    4    12 1.071022e-10 3.223569e-11
#> Dataset_13    2    13 5.073616e-16 4.282824e-16
#> Dataset_14   10    14 9.719419e-01 4.830874e-01
#> Dataset_15    6    15 9.516504e-01 1.586084e-01
#> Dataset_16   20    16 9.834545e-01 7.302771e-01
#> Dataset_17    3    17 2.183972e-11 7.280192e-12
#> Dataset_18   16    18 9.783377e-01 6.677098e-01
#> Dataset_19   12    19 9.741608e-01 5.647524e-01
#> Dataset_20    9    20 9.715404e-01 4.287703e-01
#> 
#> $message
#> [1] "5 datasets are significant at alpha level 0.10. Total datasets tested: 20."
#> 
#> $significant_units
#> [1] "Dataset_4"  "Dataset_13" "Dataset_17" "Dataset_12" "Dataset_1"
```
