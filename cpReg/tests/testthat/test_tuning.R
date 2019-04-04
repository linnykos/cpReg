context("Test tuning")

## .convert_cp_idx is correct

test_that(".convert_cp_idx works", {
  fold_id <- rep(1:5, times = 10)
  fold <- 4
  partition <- c(0, 10, length(fold_id[fold_id != 4]))
  res <- .convert_cp_idx(partition, fold_id, fold)

  expect_true(length(res) == length(partition))
  expect_true(res[1] == 0)
  expect_true(res[length(res)] == length(fold_id))
})

test_that(".convert_cp_idx for a different case", {
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)

    fold_id <- c(rep(1:5, times = 10), c(1:sample(1:5, 1)))
    fold <- sample(1:5, 1)
    partition <- c(0, sample(1:length(fold_id[fold_id != fold])-1, 1),
                   length(fold_id[fold_id != fold]))

    res <- .convert_cp_idx(partition, fold_id, fold)

    tmp <- matrix(NA, ncol = length(fold_id), nrow = 2)
    tmp[1,] <- fold_id
    tmp[2,which(fold_id != fold)] <- 1:length(fold_id[fold_id != fold])

    res2 <- c(0, which(tmp[2,] == partition[2]), which(tmp[2,] == partition[3]))

    all(res == res2)
  })

  expect_true(all(bool_vec))
})


