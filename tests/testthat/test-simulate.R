test_that("simulate_data uses a scalar sd as-is", {
  set.seed(1)
  # Regression test: sample(2, ...) shorthand used to draw sd from {1, 2}
  dat <- fashr:::simulate_data(function(x) 0 * x, x = 1:20, sd = 2)
  expect_equal(dat$sd, rep(2, 20))
  # A vector sd is sampled from the given values only
  dat2 <- fashr:::simulate_data(function(x) 0 * x, x = 1:50, sd = c(0.5, 3))
  expect_true(all(dat2$sd %in% c(0.5, 3)))
})

test_that("simulate_process returns a complete dataset for each type", {
  set.seed(1)
  for (type in c("linear", "quadratic", "nonlinear", "nondynamic")) {
    dat <- simulate_process(x = seq(0, 16, length.out = 10), type = type)
    expect_named(dat, c("x", "y", "truef", "sd"))
    expect_equal(nrow(dat), 10)
    expect_true(all(is.finite(dat$y)))
  }
})

test_that("fash_bma_sampling validates weight lengths", {
  fit <- get_test_fit()
  expect_error(
    fashr:::fash_bma_sampling(
      data_i = fit$fash_data$data_list[[1]],
      posterior_weights = c(0.5, 0.5),
      psd_values = c(0.1, 0.5, 1),
      refined_x = 1:5, M = 5,
      Si = fit$fash_data$S[[1]],
      likelihood = "gaussian", order = 2
    ),
    "same length"
  )
})
