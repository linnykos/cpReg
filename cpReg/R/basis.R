construct_AR_basis <- function(dat, thres_u = Inf, lag = 1){
  stopifnot(ncol(dat) > 0)
  if(nrow(dat) <= lag) return(matrix(0, ncol = lag*ncol(dat), nrow = 1))

  M <- ncol(dat); TT <- nrow(dat)

  res_list <- lapply(1:lag, function(x){
    rbind(matrix(0, ncol = M, nrow = x), dat[-((TT - x + 1):TT),])
  })

  transform_dat <- do.call(cbind, res_list)
  stopifnot(is.matrix(transform_dat))
  transform_dat[transform_dat > thres_u] <- thres_u
  transform_dat
}
