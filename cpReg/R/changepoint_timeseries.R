# Poisson time series, penality is L_{1,1}
cp_ar1 <- function(dat, thres_u = round(quantile(dat[dat > 0], probs = 0.75)), lambda = NA,
                   basis_function = construct_basis,
                   verbose = F, ...){

  stopifnot(is.matrix(dat), nrow(dat) > 1)
  stopifnot(all(dat >= 0), all(dat %% 1 == 0))
  stopifnot(any(dat < thres_u))

  d <- ncol(dat); n <- nrow(dat)
  transform_dat <- construct_basis(dat, thres_u = thres_u, ...)
  stopifnot(nrow(dat) == nrow(transform_dat), is.matrix(transform_dat))

  # perform row-wise glmnets
  if(verbose) print("Starting to estimate coefficients row-wise: ")
  res_list <- lapply(1:d, function(x){
    if(d > 10 && x %% floor(d/10) == 0) cat('*')
    .rowwise_glmnet(dat[-n,x], transform_dat[-n,x])
  })

  # if lambda is NA, combine the results together
  if(is.na(lambda)) {
    print("Starting to combine across different lambdas")
    est <- .combine_lambdas(res_list, dat, transform_dat, verbose = verbose)
  } else {
    est <- .extract_lambdas(res_list, lambda)
  }

  structure(list(nu = est$nu, A = est$A), class = "cp_ar1")
}

########

.objective_func <- function(){

}

.objective_func_row <- function(){

}

.nll <- function(){

}

.nll_row <- function(){

}

# see https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html#poi
.rowwise_glmnet <- function(response, covariates, lambda = NA){
  if(is.na(lambda)){
    glmnet::cv.glmnet(covariates, response, family = "poisson", alpha = 1)

  } else {
    res <- glmnet::glmnet(covariates, response, family = "poisson", alpha = 1)
    as.numeric(glmnet::coef.glmnet(res, s = lambda))
  }
}
