context("Test create data")

## create_data is correct

test_that("create_data works", {
  res <- create_data(list(c(1,1,1), c(2,-1,2)), c(0, 50, 100))

  expect_true(is.list(res))
  expect_true(all(names(res) == c("X", "y")))
  expect_true(all(dim(res$X) == c(100, 3)))
  expect_true(length(res$y) == 100)
})
