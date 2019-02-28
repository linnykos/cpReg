wbs <- function(data,
                data_length_func,
                compute_cusum_func,
                tau_function, M = 100, delta = 1,
                verbose = F, ...){
  # initialize
  n <- data_length_func(data)
  random_intervals <- .generate_intervals(n, M)
  q <- dequer::queue()
  dequer::pushback(q, c(0,n))
  b_vec <- numeric(0)
  counter <- 0

  # iterate over all intervals recursively
  while(length(q) > 0){
    if(verbose) print(counter)
    interval <- dequer::pop(q)
    if(interval[2] - interval[1] < 2*delta+1) next()
    interval_list <- .truncate(interval, random_intervals)
    res_list <- lapply(interval_list, function(x){
      .find_breakpoint(data = data, interval = x, delta = delta,
                       compute_cusum_func = compute_cusum_func,
                       verbose = verbose, ...)})
    res <- res_list[[which.max(sapply(res_list, function(x){x$val}))]]

    # if passes threshold, recurse
    if(res$val >= tau_function(data = data, interval = interval, ...)){
      b_vec <- c(b_vec, res$b)
      if(verbose) print(paste0("Pushing ", interval[1], " - ", res$b, " - ", interval[2]))
      dequer::pushback(q, c(interval[1], res$b))
      dequer::pushback(q, c(res$b+1, interval[2]))
    }

    counter <- counter+1
  }

  c(0, sort(b_vec), n)
}

################

.generate_intervals <- function(n, M){
  stopifnot(n >= 2, M >= 0, M %% 1 == 0)

  if(M > 0){
    random_intervals <- lapply(1:M, function(x){
      sort(sample(1:n, 2, replace = F))
    })
  } else {
    random_intervals <- list()
  }

  random_intervals[[M+1]] <- c(1,n)
  .remove_duplicate(random_intervals)
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

.find_breakpoint <- function(data, interval, delta = 1, compute_cusum_func,
                             verbose = F, ...){
  # check length
  stopifnot(interval[2] > interval[1])
  if(interval[2] - interval[1] < 2*delta+1) return(list(val = NA, b = NA))

  vec <- sapply((interval[1]+delta):(interval[2]-1-delta), function(x){
    compute_cusum_func(data, interval[1], interval[2], x, ...)
  })

  if(verbose) print(paste0(interval[1], " - ", which.max(vec)+interval[1]+delta-1, " - ", interval[2]))

  list(val = max(vec), b = which.max(vec)+interval[1]+delta-1)
}
