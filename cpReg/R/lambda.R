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

.combine_lambdas <- function(res_list, dat, transform_dat, verbose = verbose){
  stopifnot(all(sapply(res_list, function(x){class(x)=="cv.glmnet"})))
  stopifnot(all(sapply(res_list, function(x){"fishnet" %in% class(x$glmnet.fit)})))

  # form the sequence of lambdas first
  if(verbose) paste0("Forming sequence of lambdas")
  lambda_min <- sapply(res_list, function(x){x$lambda.1se})
  lambda_all <- sort(unique(c(lambda_min)))

  # form many estimates of A
  if(verbose) paste0("Forming and evaluating many estimates of A: ")
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

