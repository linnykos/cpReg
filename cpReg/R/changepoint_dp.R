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
                          skip_interval = 1, verbose = T){
  # create hash table
  TT <- nrow(dat); M <- ncol(dat)
  TT_seq <- seq(1, TT, by = skip_interval)
  TT_seq[length(TT_seq)] <- TT
  h <- vector("list", length(TT_seq)-1)


  for(i in TT_seq){
    if(i == 1) next()
    idxi <- which(TT_seq[-1] == i)
    if(verbose & idxi %% floor(length(TT_seq)/10) == 0) cat('*')

    obj_vec <- rep(NA, idxi-1)
    lis <- vector("list", length = length(TT_seq[TT_seq < i]))

    # search
    for(j in TT_seq[TT_seq < i]){
      idxj <- which(TT_seq == j)
      if(i - j < min_spacing){
        lis[[idxj]] <- list(A = NA)
        obj_vec[idxj] <- Inf
      } else {
        lis[[idxj]] <- stationary_ar(dat = dat[j:i, , drop = F], thres_u = thres_u,
                                  lambda = lambda * sqrt(i-j), verbose = F,
                                  intercept = F)
        tmp <- h[[idxj]]$obj_val
        if(is.null(tmp)) tmp <- 0
        obj_vec[idxj] <- lis[[idxj]]$obj_val + tmp + gamma*TT
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
    h[[idxi]] <- h_obj
  }

  h[[length(TT_seq)-1]]$partition[1] <- 0
  h[[length(TT_seq)-1]]
}


generative_model_cp <- function(nu, A_list, changepoint_perc,
                                timesteps = 10, thres_u = 5,
                                basis_function = construct_AR_basis,
                                warning_val = 1000, ...){

  stopifnot(all(changepoint_perc < 1), all(changepoint_perc > 0),
            length(changepoint_perc) == length(unique(changepoint_perc)),
            all(changepoint_perc == sort(changepoint_perc)))
  stopifnot(length(changepoint_perc)+1 == length(A_list))
  stopifnot(length(unique(sapply(A_list, ncol))) == 1,
            length(unique(sapply(A_list, nrow))) == 1)

  M <- length(nu); MK <- ncol(A_list[[1]]); TT <- timesteps
  dat <- matrix(0, ncol = M, nrow = TT)

  changepoint_idx <- round(TT*changepoint_perc)
  idx <- 1

  for(time in 1:TT){
    if(idx <= length(changepoint_idx) && time > changepoint_idx[idx]) idx <- idx+1

    obs_vec <- sapply(1:M, function(i){
      vec <- construct_AR_basis(dat[which(1:TT <= time),,drop = F], thres_u = thres_u)
      vec <- vec[nrow(vec),]
      natural_param <- as.numeric(nu[i] + A_list[[idx]][i,] %*% vec)

      if(!is.na(warning_val) & as.numeric(natural_param) > log(warning_val)){
        stop(paste0("Natural parameters are exploding, with values above log of ", warning_val))
      }
      stats::rpois(1, lambda = exp(natural_param))
    })

    dat[time,] <- obs_vec
  }

  structure(list(dat = dat, partition = changepoint_idx), class = "SEPP_cp")
}

########

.construct_list_obj <- function(obj_val, partition, A_list){
  stopifnot(length(A_list) == length(partition)-1)
  stopifnot(all(partition == sort(partition)), length(partition) == length(unique(partition)))

  list(obj_val = obj_val, partition = partition, A_list = A_list)
}

