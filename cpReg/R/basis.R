#' Construction of the basis values for auto-regressive processes
#'
#' The first few rows of the returned object will have 0's, implying that
#' the "initial values" of what is regressed on are 0.
#'
#' @param dat Count data, where each row represents a different time step and each
#' column represents a different variable
#' @param thres_u A positive threshold for the saturation effect
#' @param lag Amount of lag, should be a positive integer
#'
#' @return A matrix that has \code{nrow(dat)} rows and \code{ncol(dat)*lag}
#' columns.
#' @export
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
