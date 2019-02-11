screening <- function(fit, tau, M = 100){
  # initialize
  n <- nrow(fit)
  random_intervals <- .generate_intervals(n, M)
  q <- dequer::queue()
  dequer::pushback(q, c(1,n)) #Q: should this be 0 or 1???
  b_vec <- numeric(0)

  # iterate over all intervals recursively
  while(length(q) > 0){
    interval <- dequer::pop(q)
    interval_list <- .truncate(interval, random_intervals)
    res_list <- sapply(interval_list, .find_breakpoint, fit = fit)
    res <- res_list[[which.max(sapply(res_list, function(x){x$val}))]]

    # if passes threshold, recurse
    if(res$val >= tau){
      b_vec <- c(b_vec, res$b)
      dequer::pushback(q, c(interval[1], res$b))
      dequer::pushback(q, c(res$b+1, interval[2]))
    }
  }

  c(0, sort(b_vec), n)
}

################

.generate_intervals <- function(n, M){
  stopifnot(n >= 2)

  random_intervals <- lapply(1:M, function(x){
    sort(sample(1:n, 2, replace = F))
  })

  random_intervals[[M+1]] <- c(1,n)
  random_intervals
}

.truncate <- function(interval, random_intervals){
  stopifnot(interval[1] < interval[2])

  interval_list <- lapply(random_intervals, function(x){
    stopifnot(x[1] < x[2])
    if(x[1] >= interval[2]) return(numeric(0))
    if(interval[1] >= x[2]) return(numeric(0))

    tmp <- interval
    if(x[1] >= interval[1] & x[1] <= interval[2]) tmp[1] <- x[1]
    if(x[2] >= interval[1] & x[2] <= interval[2]) tmp[2] <- x[2]
    tmp
  })

  interval_list <- interval_list[which(sapply(interval_list, length) > 0)]
  stopifnot(length(interval_list) >= 0)
  .remove_duplicate(interval_list)
}

.remove_duplicate <- function(interval_list){
  mat <- matrix(unlist(interval_list), ncol = 2, byrow = T)
  vec <- mat[,1]*max(mat[,2]+1)+mat[,2]
  idx <- which(!duplicated(vec))

  interval_list[idx]
}

.find_breakpoint <- function(fit, interval){
  # check length
  stopifnot(interval[2] > interval[1])
  if(interval[2] - interval[1] < 3) return(NA)

  vec <- sapply((interval[1]+1):(interval[2]-1), function(x){
    .compute_cusum(fit, interval[1], interval[2], x)
  })

  list(val = max(vec), b = which.max(vec)+interval[1])
}

.compute_cusum <- function(fit, s, e, b){
  #Q: check to see if anything funny happens w/ s ane e being too short...
  stopifnot(s < b, b < e, s >= 0, e <= nrow(fit))

  factor1 <- sqrt((e-b)/((e-s)*(b-s)))
  factor2 <- sqrt((b-s)/((e-s)*(e-b)))
  .l2norm(factor1*colSums(fit[(s+1):b,,drop = F]) - factor2*colSums(fit[(b+1):e,,drop = F]))
}

.l2norm <- function(x){
  sqrt(sum(x^2))
}
