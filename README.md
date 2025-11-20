
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
#> Starting data setup...
#> Completed data setup in 0.00 seconds.
#> Starting likelihood computation...
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%  |                                                                              |=======                                                               |  10%  |                                                                              |==========                                                            |  15%  |                                                                              |==============                                                        |  20%  |                                                                              |==================                                                    |  25%  |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%  |                                                                              |============================                                          |  40%  |                                                                              |================================                                      |  45%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  55%  |                                                                              |==========================================                            |  60%  |                                                                              |==============================================                        |  65%  |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  75%  |                                                                              |========================================================              |  80%  |                                                                              |============================================================          |  85%  |                                                                              |===============================================================       |  90%  |                                                                              |==================================================================    |  95%  |                                                                              |======================================================================| 100%
#> Completed likelihood computation in 7.78 seconds.
#> Starting empirical Bayes estimation...
#> Completed empirical Bayes estimation in 0.00 seconds.
#> fash object created successfully.
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

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

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

``` r
plot(datasets[[1]]$t, datasets[[1]]$y, type = "p", col = "black", ylab = "y", xlab = "t")
lines(fitted$x, fitted$mean, col = "red")
polygon(c(fitted$x, rev(fitted$x)), c(fitted$lower, rev(fitted$upper)), col = rgb(1, 0, 0, 0.2), border = NA)
matlines(fitted$x, fitted_samps[, 1:5], col = "blue", lty = 2, lwd = 0.5)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

Compute the FDR, and highlight the datasets with FDR \< 0.1:

``` r
fdr_result <- fdr_control(fash_fit, alpha = 0.1, plot = TRUE)
#> 5 datasets are significant at alpha level 0.10. Total datasets tested: 20.
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

``` r
str(fdr_result)
#> List of 1
#>  $ fdr_results:'data.frame': 20 obs. of  2 variables:
#>   ..$ index: int [1:20] 4 13 17 12 1 15 9 3 20 14 ...
#>   ..$ FDR  : num [1:20] 3.49e-16 4.28e-16 7.28e-12 3.22e-11 3.99e-09 ...
```
