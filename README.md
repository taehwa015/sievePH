# `smlePH`: R package for fitting the proportional hazards model with sieve maximum full likelihood estimation

## Overview

`smlePH` is a R package to fit the proportional hazards model with the sieve maximum full likelihood estimation. Unlike the conventional partial likelihood approach, the full likelihood method utilizes whole likelihood information in estimating procedure. In this package, we provide two functions `smle_ph` and `smle_resid`. `smle_ph` fits the proportional hazards model and estimates the baseline cumulative hazard function. `smle_resid` estimates two types of residuals, (i) score residual for each covariate and (ii) deviance residual for linear combination of all covariates.

## Installation


The released version of `smlePH` package from CRAN is available in CRAN
```r
install.packages('smlePH')
```
and the development version is also available in Github 
```r
devtools::install_github(repo='taehwa015/smlePH')
```



