#' Poisson time series, penality is L_{1,1}
#'
#' Solves the optimization problem (2.8) on https://arxiv.org/pdf/1802.04838.pdf.
#' The penalty is fixed to be the L_{1,1} (sum of absolute value), but the
#' construction of the basis can be passed in as a function.
#'
#' The implementation of the function fits \code{ncol(dat)} Lasso fits under
#' the Poisson model and concatenates the solutions together.
#' If lambda is not set, then a heurstical way to concetenate solutions via
#' cross-validation is used.
#'
#' @param dat Count data, where each row represents a different time step and each
#' column represents a different variable
#' @param thres_u A positive threshold for the saturation effect
#' @param lambda Tuning parameter for the regression problem, allowed to be \code{NA}
#' @param basis_function The function to generate the variables to regress onto, which
#' can possibly take addition inputs
#' @param intercept Boolean to allow for estimating intercepts
#' @param verbose Boolean
#' @param ... Addition parameters for the basis function
#'
#' @return A list containing the final-used \code{lambda} value, the final
#' objective value (\code{obj_val}), the estimated intercepts (\code{nu}) and
#' the estimated transition matrix (\code{A}).
#' @export
stationary_ar <- function(dat, thres_u = round(stats::quantile(dat[dat > 0], probs = 0.75)),
                          lambda = NA, sparsity = NA,
                          basis_function = construct_AR_basis, intercept = T,
                          verbose = F, ...){

  stopifnot(is.matrix(dat), nrow(dat) > 1)
  stopifnot(all(dat >= 0), all(dat %% 1 == 0))
  stopifnot(any(dat < thres_u), thres_u >= 0)

  M <- ncol(dat); TT <- nrow(dat)
  transform_dat <- construct_AR_basis(dat, thres_u = thres_u, ...)
  stopifnot(nrow(dat) == nrow(transform_dat), is.matrix(transform_dat))

  # perform row-wise glmnets
  if(verbose) cat("\nStarting to estimate coefficients row-wise: ")
  res_list <- lapply(1:M, function(x){
    if(M > 10 && x %% floor(M/10) == 0) cat('*')
    .rowwise_glmnet(dat[,x], transform_dat, (is.na(lambda) & is.na(sparsity)), intercept =  intercept)
  })

  # if lambda is NA, combine the results together
  if(verbose) cat("\nStarting to assemble estimator")
  if(is.na(lambda)) {
    if(!is.na(sparsity)){
      est <- .extract_lambdas_sparsity(res_list, sparsity, verbose = verbose)
    } else {
      est <- .combine_lambdas(res_list, dat, transform_dat, verbose = verbose)
    }
  } else {
    est <- .extract_lambdas(res_list, lambda)
  }

  obj_val <- .objective_func(est$nu, est$A, dat, transform_dat, lambda = est$lambda)
  stopifnot(nrow(est$A) == ncol(dat), ncol(est$A) == ncol(transform_dat),
            length(est$nu) == ncol(dat))
  structure(list(lambda = est$lambda, obj_val = obj_val, nu = est$nu, A = est$A), class = "cp_ar1")
}

#' Generate a dataset from a true model
#'
#' Samples from the model in (2.1), (2.2) with saturation in (2.6) in
#' https://arxiv.org/pdf/1802.04838.pdf. The
#' construction of the basis can be passed in as a function.
#'
#' @param nu Vector in intercepts, one for each variable
#' @param A Square transition matrix
#' @param timesteps Number of time steps to proceed forward with
#' @param thres_u A positive threshold for the saturation effect
#' @param basis_function The function to generate the variables to regress onto, which
#' can possibly take addition inputs
#' @param warning_val If not \code{NA}, the code will stop if any of the
#' parameters of the Poisson exceed this number. The default is set to 1000.
#' @param ... Addition parameters for the basis function
#'
#' @return A matrix that has \code{timesteps} rows and \code{ncol(A)} columns
#' of non-negative integers
#' @export
#'
generative_model <- function(nu, A, timesteps = 10, thres_u = 5,
                             basis_function = construct_AR_basis,
                             warning_val = 1000, ...){
  stopifnot(nrow(A) == length(nu))

  M <- length(nu); MK <- ncol(A); TT <- timesteps
  dat <- matrix(0, ncol = M, nrow = TT)

  for(time in 1:TT){
    obs_vec <- sapply(1:M, function(i){
      vec <- construct_AR_basis(dat[which(1:TT < time),,drop = F], thres_u = thres_u, ...)
      vec <- vec[nrow(vec),]
      natural_param <- as.numeric(nu[i] + A[i,] %*% vec)

      if(!is.na(warning_val) & as.numeric(natural_param) > log(warning_val)){
        stop(paste0("Natural parameters are exploding, with values above log of ", warning_val))
      }
      stats::rpois(1, lambda = exp(natural_param))
    })

    dat[time,] <- obs_vec
  }

  dat
}

########

.objective_func <- function(nu, A, dat, transform_dat, lambda){
  .nll(nu, A, dat, transform_dat) + lambda *sum(abs(A))
}

.objective_func_row <- function(nu_val, A_vec, dat_vec, transform_dat, lambda){
  .nll_row(nu_val, A_vec, dat_vec, transform_dat) + lambda * sum(abs(A_vec))
}

.nll <- function(nu, A, dat, transform_dat){
  stopifnot(nrow(dat) == nrow(transform_dat), nrow(A) == ncol(dat),
            ncol(A) == ncol(transform_dat))
  TT <- nrow(dat)
  intercept_mat <- t(sapply(1:(TT-1), function(x){nu}))

  sum(exp(intercept_mat+ transform_dat[1:(TT-1),]%*%t(A))) - sum(dat[2:TT,] *(intercept_mat + transform_dat[1:(TT-1),]%*%t(A)))
}

# for each dimension
.nll_row <- function(nu_val, A_vec, dat_vec, transform_dat){
  stopifnot(length(dat_vec) == nrow(transform_dat))
  TT <- length(dat_vec)

  sum(exp(nu_val + A_vec %*% t(transform_dat[1:(TT-1),]))) - sum(dat_vec[2:TT] * (nu_val + A_vec %*% t(transform_dat[1:(TT-1),])))
}

# see https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html#poi
.rowwise_glmnet <- function(response, covariates, cv = TRUE, intercept = TRUE){
  stopifnot(length(response) == nrow(covariates))
  TT <- length(response)
  stopifnot(TT > 2)

  if(cv){
    if(TT < 30) stop("Not enough observations to perform cross validation")
    glmnet::cv.glmnet(x = covariates[1:(TT-1),,drop = F], y = response[2:TT], family = "poisson", alpha = 1,
                      nfolds = min(max(floor(TT/10), 3), 10), intercept = intercept)

  } else {
    glmnet::glmnet(x = covariates[1:(TT-1),,drop = F], y = response[2:TT], family = "poisson", alpha = 1,
                   intercept = intercept)
  }
}
