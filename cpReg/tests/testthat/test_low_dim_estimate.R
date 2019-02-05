context("Test low dimensional estimator")

## .format_sequence is correct

test_that(".format_sequence works", {
  res <- .format_sequence(10, 2)
  expect_true(all(res == c(0,2,4,6,8,10)))
})

test_that(".format_sequence works if n isn't evenly divided by delta", {
  res <- .format_sequence(13, 5)
  expect_true(all(res == c(0,5,13)))
})

#######

## .construct_list_obj is correct

test_that(".construct_list_obj works", {
  res <- .construct_list_obj(10, c(1,10), 15)
  expect_true(is.list(res))
  expect_true(all(names(res) == c("obj_val", "partition", "coef_list")))
})

########

## .local_regression is correct

test_that(".local_regression works", {
  set.seed(10)
  X <- matrix(rnorm(200), nrow = 20)
  y <- rnorm(20)
  res <- .local_regression(X,y)

  expect_true(is.list(res))
  expect_true(all(names(res) == c("obj_val", "coef")))
})

#############

## .special_concatenate is correct

test_that(".special_concatenate works with is_list = F", {
  res <- .special_concatenate(1:5, 11:15, F)
  expect_true(all(res == c(1:5,11:15)))
})

test_that(".special_concatenate works with is_list = T", {
  res <- .special_concatenate(list(1:5, 11:15), 21:25, T)
  expect_true(is.list(res))
  expect_true(length(res) == 3)
})

