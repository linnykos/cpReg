#' Changepoint multivariate Poisson time series, penality is L_{1,1}
#'
#' Does a screening based on the multidimensional CUSUM.
#'
#' @param fit The output object from \code{changepoint_dp()}
#' @param gamma Tuning parameter to control the number of output changepoints. Cannot be \code{NA}
#'
#' @return a list with a vector \code{partition} denoting the resulting
#' partition and a vector \code{cusum_vec} denoting the resulting CUSUM at each of the middle indices
#' in \code{partition}
#' @export
changepoint_screening <- function(fit, gamma){
  partition <- fit$partition
  A_list <- .unravel(partition, fit$A_list)

  .enumerate_all_cp <- function(partition){
    len <- length(partition)
    if(len < 3) return(NA)
    for(i in 1:(len-2)){
      cusum <- .multidimen_cusum(A_list[(partition[i]+1):partition[i+1]],
                                 A_list[(partition[i+1]+1):partition[i+2]])
      val <- .l2norm(cusum)
      if(val < gamma) return(i+1)
    }

    NA
  }

  while(TRUE){
    tmp <- .enumerate_all_cp(partition)
    if(is.na(tmp)) break()
    partition <- partition[-tmp]
  }

  # output the cusum vector
  if(length(partition) > 2){
    cusum_vec <- sapply(1:(length(partition)-2), function(i){
      .l2norm(.multidimen_cusum(A_list[(partition[i]+1):partition[i+1]],
                                A_list[(partition[i+1]+1):partition[i+2]]))
    })
  } else {
    cusum_vec <- NA
  }


  list(partition = partition, cusum_vec = cusum_vec)
}

# converts k-1 estimates of A into T estimates of A
.unravel <- function(partition, A_list){
  stopifnot(partition[1] == 0, all(sort(partition) == partition),
            all(partition == unique(partition)),
            length(partition) == length(A_list)+1)

  len <- length(partition)
  TT <- partition[len]
  lis <- vector("list", TT)
  counter <- 1
  for(i in 1:TT){
    if(i > partition[counter+1]) counter <- counter+1
    lis[[i]] <- A_list[[counter]]
  }

  lis
}

.multidimen_cusum <- function(list1, list2){
  len1 <- length(list1)
  len2 <- length(list2)

  sqrt(len2/((len1+len2)*len1))*Reduce('+', list1) -
    sqrt(len1/((len1+len2)*len2))*Reduce('+', list2)
}

.l2norm <- function(x){
  sqrt(sum(as.numeric(x)^2))
}
