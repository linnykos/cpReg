rm(list=ls())
M <- 100
tau <- 10
fit <- matrix(0, nrow = 200, ncol = 3)
fit[1:100,] <- c(rep(1,100), rep(0, 100), rep(1,100))
fit[101:200,] <- c(rep(-1,100), rep(0, 100), rep(-1,100))

set.seed(10)
n <- nrow(fit)
random_intervals <- .generate_intervals(n, M)
q <- dequer::queue()
dequer::pushback(q, c(0,n)) #Q: should this be 0 or 1???
b_vec <- numeric(0)

counter <- 0

# iterate over all intervals recursively
while(length(q) > 0){
  print(counter)
  interval <- dequer::pop(q)
  interval_list <- .truncate(interval, random_intervals)
  res_list <- lapply(interval_list, .find_breakpoint, fit = fit)
  res <- res_list[[which.max(sapply(res_list, function(x){x$val}))]]

  # if passes threshold, recurse
  if(res$val >= tau){
    b_vec <- c(b_vec, res$b)
    dequer::pushback(q, c(interval[1], res$b))
    dequer::pushback(q, c(res$b+1, interval[2]))
  }

  counter <- counter+1
}
