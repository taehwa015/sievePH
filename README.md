# `sievePH`: R package for fitting the proportional hazards model with sieve maximum full likelihood estimation


## Overview

`sievePH` is the R package to fit the proportional hazards model with sieve maximum full likelihood estimation.

## Installation
```r
devtools::install_github(repo='taehwa015/sievePH')
```

## Status of development

The code provided here is only being provided for research purposes.

## Documentation

Vignette is available at [here](http://htmlpreview.github.io/?https://github.com/taehwa015/sievePH/blob/master/vignettes/sievePH.html).

## Usage

The `sievePH` package provides a sieve maximum full likelihood estimation for the proportional hazards model.
See Halabi et al. (2024+) for more detailed description of the method.


Below example is the phase 3 metastatic colorectal cancer clinical trial study.
Since the event times are possibly correlated within same patient, 
we regard patient id as the cluster. 
To adjust informative cluster sizes, we further consider weight function,
motivated by Wang and Zhao (2008).
Larger cluster will be underweighted by letting `alpha = 1`, 
while cluster structure will be ignored `alpha = 0`.
Note that by letting `id = NULL`, we fit the univariate AFT model.
```r
library(PICBayes)
data(mCRC)
dt0 = as.data.frame(mCRC)
d = with(dt0,
         data.frame(U = ifelse(is.na(L), 0, L),
                    V = ifelse(is.na(R), Inf, R),
                    Delta = 1-IC,
                    x1 = TRT_C,
                    x2 = KRAS_C,
                    id = SITE))
U = d$U; V = d$V; X = cbind(d$x1, d$x2); Delta = d$Delta; id = d$id
aft_rank(U = U, V = V, X = X, Delta = Delta, id = id, 
         alpha = 1, type = "gehan", R = 10)
aft_rank(U = U, V = V, X = X, Delta = Delta, id = id, 
         alpha = 1, type = "logrank", R = 10)
```

## Reference


