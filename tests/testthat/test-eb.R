test_that("fash_eb_est returns normalized prior and posterior weights", {
  set.seed(1)
  grid <- c(0, 0.5, 1)
  # Half the datasets strongly prefer the null, half the last component
  L <- rbind(
    matrix(rep(c(0, -20, -30), each = 10), nrow = 10),
    matrix(rep(c(-30, -20, 0), each = 10), nrow = 10)
  )
  result <- fash_eb_est(L, penalty = 1, grid = grid)
  expect_equal(sum(result$prior_weight$prior_weight), 1, tolerance = 1e-6)
  expect_true(all(abs(rowSums(result$posterior_weight) - 1) < 1e-8))
  w_null <- result$prior_weight$prior_weight[result$prior_weight$psd == 0]
  expect_equal(w_null, 0.5, tolerance = 0.05)
})

test_that("fash_eb_est with penalty > 1 is robust to strongly negative log-likelihoods", {
  set.seed(1)
  grid <- c(0, 0.5, 1)
  # Log-likelihoods around -5000: exp() would underflow without rescaling
  L <- rbind(
    matrix(rep(c(-5000, -5020, -5030), each = 10), nrow = 10),
    matrix(rep(c(-5030, -5020, -5000), each = 10), nrow = 10)
  )
  result <- fash_eb_est(L, penalty = 10, grid = grid)
  expect_true(all(is.finite(result$prior_weight$prior_weight)))
  expect_true(all(is.finite(result$posterior_weight)))
  expect_true(all(abs(rowSums(result$posterior_weight) - 1) < 1e-8))
  # The Dirichlet penalty pulls extra weight onto the null component
  result_nopen <- fash_eb_est(L, penalty = 1, grid = grid)
  w_null_pen <- result$prior_weight$prior_weight[result$prior_weight$psd == 0]
  w_null_nopen <- result_nopen$prior_weight$prior_weight[result_nopen$prior_weight$psd == 0]
  expect_gte(w_null_pen, w_null_nopen - 1e-8)
})

test_that("fash_post_ordering orders by the requested metric", {
  set.seed(1)
  grid <- c(0, 0.5, 1)
  L <- matrix(rnorm(15), nrow = 5)
  eb <- fash_eb_est(L, penalty = 1, grid = grid)
  res <- fash_post_ordering(eb, ordering = "mean")
  psd_values <- as.numeric(colnames(eb$posterior_weight))
  means <- as.numeric(eb$posterior_weight %*% psd_values)
  expect_equal(res$ordered_indices, order(means))
  expect_equal(res$ordered_metrics, sort(means))
})
