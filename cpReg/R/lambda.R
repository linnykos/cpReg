.lambda_oracle <- function(obj,
                          thres_u = round(stats::quantile(obj$dat[obj$dat > 0], probs = 0.75)),
                          basis_function = construct_AR_basis, intercept = T){

  stopifnot(class(obj) == "SEPP_cp")
  changepoint_idx <- c(0, obj$partition, nrow(obj$dat))
  k <- length(changepoint_idx)-1
  lambda_vec <- sapply(1:k, function(i){
    len <- changepoint_idx[i+1]-changepoint_idx[i]

    lambda <- stationary_ar(obj$dat[(changepoint_idx[i]+1):changepoint_idx[i+1],],
                            thres_u = thres_u, lambda = NA,
                            basis_function = basis_function, intercept = intercept,
                            verbose = F)$lambda
    lambda/sqrt(len) * (1/len)
  })

  sum(lambda_vec)
}

#####

.extract_lambdas <- function(res_list, lambda){
  stopifnot(all(sapply(res_list, function(x){"fishnet" %in% class(x)})))

  res <- lapply(res_list, function(x){
    as.numeric(glmnet::coef.glmnet(x, s = lambda))
  })

  stopifnot(length(unique(sapply(res, length))) == 1)
  res <- do.call(rbind, res)

  nu <- res[,1]
  A <- res[,-1]

  list(lambda = lambda, nu = nu, A = A)
}

.extract_lambdas_sparsity <- function(res_list, s, verbose = F, max_length = 100){
  stopifnot(all(sapply(res_list, function(x){"fishnet" %in% class(x)})))
  stopifnot(0 <= s, s <= 1)

  lambda_all <- sort(unique(unlist(lapply(res_list, function(x){x$lambda}))))
  if(length(lambda_all) > max_length){
    lambda_all <- lambda_all[sort(unique(round(seq(1, length(lambda_all), length.out = max_length))))]
  }

  if(verbose) cat(paste0("\nEvaluating many estimates of A for sparsity: "))
  sparsity_vec <- sapply(1:length(lambda_all), function(x){
    if(verbose & x %% floor(length(lambda_all)/10) == 0) cat('*')
    res <- .extract_lambdas(res_list, lambda = lambda_all[x])
    length(which(res$A == 0))/prod(dim(res$A))
  })

  lambda_chosen <- lambda_all[which.min(abs(sparsity_vec - s))]

  .extract_lambdas(res_list, lambda = lambda_chosen)
}

.combine_lambdas <- function(res_list, dat, transform_dat, verbose = verbose){
  stopifnot(all(sapply(res_list, function(x){class(x)=="cv.glmnet"})))
  stopifnot(all(sapply(res_list, function(x){"fishnet" %in% class(x$glmnet.fit)})))

  # form the sequence of lambdas first
  if(verbose) cat(paste0("\nForming sequence of lambdas"))
  lambda_min <- sapply(res_list, function(x){x$lambda.1se})
  lambda_all <- sort(unique(c(lambda_min)))

  # form many estimates of A
  if(verbose) cat(paste0("\nEvaluating many estimates of A for cv: "))
  d <- length(res_list)
  len <- length(lambda_all)
  obj_min <- Inf
  lambda_best <- numeric(0)
  A_best <- numeric(0)

  for(i in 1:length(lambda_all)){
    if(i > 10 && i %% floor(len/10) == 0) cat('*')

    A <- t(sapply(1:d, function(y){
      as.numeric(glmnet::coef.cv.glmnet(res_list[[y]], s = lambda_all[i]))
    }))

    obj <- .nll(A[,1], A[,-1], dat, transform_dat)
    if(obj < obj_min) {obj_min <- obj; A_best <- A; lambda_best <- lambda_all[i]}
  }

  list(lambda = lambda_best, nu = A_best[,1], A = A_best[,-1])
}

