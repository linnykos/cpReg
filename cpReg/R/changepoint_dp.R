# ignore ALL changepoints of size = 2
.changepoint_dp <- function(dat, thres_u = round(stats::quantile(dat[dat > 0], probs = 0.75)),
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
        tmp <- h[[j]]$obj_val #SOMETHING WERID
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
}

.construct_list_obj <- function(obj_val, partition, A_list){
  stopifnot(length(A_list) == length(partition)-1)
  stopifnot(all(partition == sort(partition)), length(partition) == length(unique(partition)))

  list(obj_val = obj_val, partition = partition, A_list = A_list)
}
