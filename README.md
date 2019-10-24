
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sclr

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/sclr)](https://cran.r-project.org/package=sclr)
[![Travis build
status](https://travis-ci.org/khvorov45/sclr.svg?branch=master)](https://travis-ci.org/khvorov45/sclr)
[![codecov](https://codecov.io/gh/khvorov45/sclr/branch/master/graph/badge.svg)](https://codecov.io/gh/khvorov45/sclr)
<!-- badges: end -->

The goal of sclr is to fit the scaled logit model from Dunning (2006)
using the maximum likelihood method.

## Installation

Install the [CRAN](https://CRAN.R-project.org) version with

``` r
install.packages("sclr")
```

Or the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("khvorov45/sclr")
```

## Model

For log likelihood, scores and second derivatives see `vignette("Math",
"sclr")`. Documentation of the main fitting function `?sclr` has details
on how the model is fit.

## Example

Usage is similar to other model fitting functions like `lm`.

``` r
library(sclr)
ideal_data <- sclr_ideal_data(seed = 20191021)
fit <- sclr(status ~ logHI, ideal_data)
summary(fit)
#> Call: status ~ logHI
#> 
#> Parameter estimates
#>     lambda     beta_0 beta_logHI 
#>  0.4633844 -5.1399247  2.0679708 
#> 
#> 95% confidence intervals
#>                 2.5 %    97.5 %
#> lambda      0.4078069  0.518962
#> beta_0     -7.2061683 -3.073681
#> beta_logHI  1.3758545  2.760087
```

For more details see `vignette("Usage", "sclr")`.

## References

Dunning AJ (2006). “A model for immunological correlates of protection.”
Statistics in Medicine, 25(9), 1485-1497. doi: 10.1002/sim.2282.
