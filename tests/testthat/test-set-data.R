test_that("fash_set_data handles matrix input", {
  set.seed(1)
  Y <- matrix(rnorm(20), nrow = 4, ncol = 5)
  smooth_var <- matrix(runif(20), nrow = 4, ncol = 5)
  S <- c(0.5, 0.8, 1.2, 1.0, 0.9)
  Omega <- diag(5)
  result <- fash_set_data(Y = Y, smooth_var = smooth_var, offset = 1, S = S, Omega = Omega)
  expect_length(result$data_list, 4)
  expect_equal(result$data_list[[2]]$y, Y[2, ])
  expect_equal(result$data_list[[2]]$x, smooth_var[2, ])
  expect_equal(result$data_list[[2]]$offset, rep(1, 5))
  expect_length(result$S, 4)
  expect_equal(as.numeric(result$S[[3]]), S)
  expect_length(result$Omega, 4)
  expect_equal(result$Omega[[4]], Omega)
})

test_that("fash_set_data handles data_list input with column names", {
  data_list <- list(
    data.frame(y = rnorm(5), x = 1:5, off = 2, sd = 0.5),
    data.frame(y = rnorm(5), x = 1:5, off = 2, sd = 0.8)
  )
  result <- fash_set_data(data_list = data_list, Y = "y", smooth_var = "x",
                          offset = "off", S = "sd")
  expect_length(result$data_list, 2)
  expect_equal(result$data_list[[1]]$offset, rep(2, 5))
  expect_equal(result$S[[2]], rep(0.8, 5))
})

test_that("fash_set_data validates offset in data_list mode", {
  data_list <- list(data.frame(y = rnorm(5), x = 1:5))
  # scalar offset is fine
  result <- fash_set_data(data_list = data_list, Y = "y", smooth_var = "x", offset = 0)
  expect_equal(result$data_list[[1]]$offset, rep(0, 5))
  # missing offset column gives an informative error
  expect_error(
    fash_set_data(data_list = data_list, Y = "y", smooth_var = "x", offset = "nope"),
    "not found"
  )
  # a non-scalar numeric offset is rejected
  expect_error(
    fash_set_data(data_list = data_list, Y = "y", smooth_var = "x", offset = c(1, 2)),
    "column name or a scalar"
  )
})

test_that("fash_set_tmbdat expands scalar S and validates dimensions", {
  data_i <- data.frame(y = rnorm(6), x = seq(0, 1, length.out = 6), offset = 0)
  tmbdat <- fash_set_tmbdat(data_i, Si = 0.5, num_basis = 5, order = 2)
  expect_equal(tmbdat$S, rep(0.5, 6))
  expect_equal(ncol(tmbdat$X), 2)
  expect_equal(ncol(tmbdat$B), 4)
  expect_error(fash_set_tmbdat(data_i, Si = c(1, 2), num_basis = 5),
               "same length")
  expect_error(fash_set_tmbdat(data_i, Omegai = diag(3), num_basis = 5),
               "same dimensions")
})
