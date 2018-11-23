#' Changepoint multivariate Poisson time series, penality is L_{1,1}
#'
#' This method ignores ALL trying changepoints of size less than \code{min_spacing}.
#' The returned element in the list, \code{partition} denotes the LAST
#' index of each interval.
#'
#' @param dat Count data, where each row represents a different time step and each
#' column represents a different variable
#' @param thres_u A positive threshold for the saturation effect
#' @param lambda Tuning parameter for the regression problem. Cannot be \code{NA}
#' @param gamma Tuning parameter to control the number of output changepoints. Cannot be \code{NA}
#' @param min_spacing Positive integer to dictate the smallest gap in changepoints. Needed for numerical stability
#' when calling \code{glmnet}
#' @param verbose Boolean
#'
#' @return A list containing \code{obj_val} (the objective value),
#' \code{partition} (the estimated partition) and \code{A_list} the estimated
#' adjecency matrices within each partition
#' @export
changepoint_dp <- function(dat, thres_u = round(stats::quantile(dat[dat > 0], probs = 0.75)),
                            lambda, gamma, min_spacing = 10,
                            verbose = T){
  # create hash table
  TT <- nrow(dat); M <- ncol(dat)
  h <- vector("list", TT)

  for(i in 2:TT){
    if(verbose & i %% floor(TT/10) == 0) cat('*')

    obj_vec <- rep(NA, i-1)
    lis <- vector("list", length = i-1)

    # search
    for(j in 1:(i-1)){
      if(i - j <= min_spacing){
        lis[[j]] <- list(A = NA)
        obj_vec[j] <- Inf
      } else {
        lis[[j]] <- stationary_ar(dat = dat[j:i, , drop = F], thres_u = thres_u,
                                  lambda = lambda * sqrt(i-j), verbose = F,
                                  intercept = F)
        tmp <- h[[j]]$obj_val
        if(is.null(tmp)) tmp <- 0
        obj_vec[j] <- lis[[j]]$obj_val + tmp + gamma
      }
    }

    # select
    idx <- which.min(obj_vec)

    # finalize
    h_prev <-  h[[idx]]
    if(all(is.null(h_prev))){
      partition <- c(1, i); A_list <- list(lis[[idx]]$A)
    } else {
      partition <- c(h_prev$partition, i)
      A_list <- h_prev$A_list; A_list[[length(A_list)+1]] <- lis[[idx]]$A
    }

    h_obj <- .construct_list_obj(min(obj_vec), partition = partition,
                                   A_list = A_list)
    h[[i]] <- h_obj
  }

  h[[TT]]$partition[1] <- 0
  h[[TT]]
}

.construct_list_obj <- function(obj_val, partition, A_list){
  stopifnot(length(A_list) == length(partition)-1)
  stopifnot(all(partition == sort(partition)), length(partition) == length(unique(partition)))

  list(obj_val = obj_val, partition = partition, A_list = A_list)
}
