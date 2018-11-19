.changepoint_screening <- function(dat, fit, gamma){
  h <- hash::hash()
  partition <- fit$partition
  A_list <- .unravel(partition, fit$A_list)

  .enumerate_all_cp <- function(partition){
    len <- length(partition)
    if(len < 3) return(NA)
    for(i in 1:(len-2)){
      cusum <- .multidimen_cusum(A_list[(partition[i]+1):partition[i+1]],
                                 A_list[(partition[i+1]+1):partition[i+2]])
      val <- .l2norm(cusum)
      if(val < gamma){
        return(i+1)
      }
    }

    NA
  }

  while(TRUE){
    tmp <- .enumerate_all_cp(partition)
    if(is.na(tmp)) break()
    partition <- partition[-tmp]
  }

  partition
}



.construct_hash_name <- function(s,t,e){
  stopifnot(s < t, t < e)

  paste0(s, "-", t, "-", e)
}

# converts k-1 estimates of A into T estimates of A
.unravel <- function(partition, A_list){
  stopifnot(partition[1] == 1, all(sort(partition) == partition),
            all(partition == unique(partition)),
            length(partition) == length(A_list)+1)

  len <- length(partition)
  TT <- partition[len]
  lis <- vector("list", len)
  counter <- 1
  for(i in 1:len){
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
