.find_breakpoint <- function(data, interval, delta = 1,
                             max_candidates = 50,
                             data_length_func,
                             compute_cusum_func,
                             verbose = F, ...){
  # check length
  stopifnot(interval[2] > interval[1])
  if(interval[2] - interval[1] < 2*delta+1) return(list(val = NA, b = NA))
  n <- data_length_func(data)

  seq_vec <- (interval[1]+delta):(interval[2]-1-delta)

  # if too many candidates, cut in half
  if(!is.na(max_candidates) && length(seq_vec) > max_candidates){
    seq_vec <- sort(intersect(seq_vec, seq(delta, n - delta, by = delta)))
  }

  if(length(seq_vec) == 0){
    return(list(val = NA, b = NA))
  }

  vec <- sapply(seq_vec, function(x){
    compute_cusum_func(data, interval[1], interval[2], x, ...)
  })

  if(verbose) print(paste0(interval[1], " - ", which.max(vec)+interval[1]+delta-1, " - ", interval[2]))

  list(val = max(vec), b = which.max(vec)+interval[1]+delta-1)
}

########

.glmnet_to_cp <- function(lambda, len){
  lambda*sqrt(len)
}

.cp_to_glmnet <- function(lambda, len){
  lambda/sqrt(len)
}
