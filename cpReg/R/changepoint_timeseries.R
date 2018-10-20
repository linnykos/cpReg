# Poisson time series, penality is L_{1,1}
cp_ar1 <- function(dat, thres_u, lambda = NA,
                   basis_function = construct_basis, ...){

  stopifnot(is.matrix(dat), nrow(dat) > 1)
  stopifnot(all(dat >= 0), all(dat %% 1 == 0))
  stopifnot(any(dat < thres_u))

  d <- ncol(dat); n <- nrow(dat)
  transform_dat <- construct_basis(dat, thres_u = thres_u, ...)
  stopifnot(nrow(dat) == nrow(transform_dat), is.matrix(transform_dat))

  # perform row-wise glmnets
  res_list <- lapply(2:n, function(x){
    .rowwise_glmnet(dat[x,], transform_dat[x,])
  })

  # if lambda is NA, combine the results together
  if(is.na(lambda)) {
    est <- .combine_lambdas(res_list)
  } else {
    est <- .extract_lambdas(res_list)
  }

  structure(list(nu = est$nu, A = est$A), class = "cp_ar1")
}

########

.objective_func <- function(){

}

.objective_func_row <- function(){

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
