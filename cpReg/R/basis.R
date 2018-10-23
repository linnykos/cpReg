construct_AR_basis <- function(dat, thres_u = Inf, lag = 1){
  M <- ncol(dat); TT <- nrow(dat)

  res_list <- lapply(1:lag, function(x){
    rbind(matrix(0, ncol = M, nrow = x), dat[-((TT - x + 1):TT),])
  })

  transform_dat <- do.call(cbind, res_list)
  transform_dat[transform_dat > thres_u] <- thres_u
  transform_dat
}
