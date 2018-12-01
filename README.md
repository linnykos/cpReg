Installation
-------
``` r
devtools::install_github("linnylin92/cpReg", subdir = "cpReg")
```

Example
-------

``` r
rm(list=ls())
library(cpReg)

## simple example to demonstrate how to set up the algorithm
# set up the data
set.seed(10)
M <- 5; TT <- 1000
nu <- rep(0.1, M)
A <- 0.3*diag(M)
obj <- cpReg::generative_model(nu, A, TT, lag = 1, thres_u = 5)

# fit the model
fit <- cpReg::stationary_ar(obj$dat, thres_u = 5, basis_function = construct_AR_basis,
                            lambda = NA, verbose = T, lag = 1)

# compute the Forbenius difference between fitted A and true A
sum((A - fit$A)^2)
``` 
