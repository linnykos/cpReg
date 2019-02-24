context("Test high dimensional feasible estimator")

## .lasso_regression is correct

test_that(".lasso_regression works", {
  set.seed(10)
  dat <- create_data(list(c(1,1,1), c(2,-1,2)), c(0, 50, 100))
  res <- .lasso_regression(dat$X, dat$y, lambda = 1)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 3)
})

##########

## .regression_cusum is correct

test_that(".regression_cusum works", {
  set.seed(10)
  dat <- create_data(list(c(1,1,1), c(2,-1,2)), c(0, 50, 100))
  res <- .regression_cusum(dat$X, dat$y, start = 1,
                           end = 100, lambda = 1, breakpoint = 50)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 3)
})

test_that(".regression_cusum is maximized appropriately", {
  set.seed(10)
  dat1 <- create_data(list(c(-10,-10,-10), c(20,20,20)), c(0, 50, 100))
  dat2 <- create_data(list(c(-10,-10,-10), c(20,20,20)), c(0, 75, 100))
  res1 <- sapply(5:95, function(x){
    .l2norm(.regression_cusum(dat1$X, dat1$y, start = 1,
                      end = 100, lambda = 1, breakpoint = x))
  })
  res2 <- sapply(5:95, function(x){
    .l2norm(.regression_cusum(dat2$X, dat2$y, start = 1,
                              end = 100, lambda = 1, breakpoint = x))
  })

  expect_true(abs(which.max(res1) - 45) < abs(which.max(res2) - 45))
  expect_true(abs(which.max(res2) - 70) < abs(which.max(res1) - 70))
})

############

## .create_node is correct

test_that(".create_node makes a valid node", {
  res <- .create_node(2, 5)
  expect_true(is_valid(res))
})

test_that(".create_node errors if start is larger than end", {
  res <- expect_error(.create_node(5,2))
})

test_that(".create_node creates the right name", {
  res <- .create_node(2, 5)
  expect_true(res$name == "2-5")
})

#############

## .get_leaves_names is correct

test_that(".get_leaves_names works", {
  data(acme, package = "data.tree")
  res <- .get_leaves_names(acme)
  expect_true(all(res == c("Go agile", "New Accounting Standards", "New Labs",
                           "New Product Line", "New Software", "Outsource", "Switch to R")))
})

###############################

## .find_leadingBreakpoint is correct

test_that(".find_leadingBreakpoint will work if there is only one leaf", {
  tree <- .create_node(1, 10)
  tree$cusum <- 10

  res <- .find_leadingBreakpoint(tree)
  expect_true(res == "1-10")
})

#####################################

## .split_node is correct

test_that(".split_node correctly returns left and right", {
  tree <- .create_node(1, 10, 5)
  res <- .split_node(tree)

  expect_true(length(res) == 2)
  expect_true(res$left$start == 1)
  expect_true(res$left$end == 5)
  expect_true(res$right$start == 6)
  expect_true(res$right$end == 10)
})

