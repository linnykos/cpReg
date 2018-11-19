.changepoint_screening <- function(dat, fit, gamma){

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
