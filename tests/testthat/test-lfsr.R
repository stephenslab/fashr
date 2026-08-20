test_that("compute_posterior_sign_prob matches the normal CDF", {
  mu <- c(-1, 0, 2)
  sigma2 <- c(1, 4, 0.25)
  res <- fashr:::compute_posterior_sign_prob(mu, sigma2)
  expect_equal(res$pos_prob, 1 - pnorm(-mu / sqrt(sigma2)))
  expect_equal(res$neg_prob, pnorm(-mu / sqrt(sigma2)))
  expect_equal(res$lfsr, pmin(res$pos_prob, res$neg_prob))
})

test_that("compute_posterior_sign_prob handles point masses (sigma = 0)", {
  res <- fashr:::compute_posterior_sign_prob(c(-1, 0, 2), c(0, 0, 0))
  expect_equal(res$pos_prob, c(0, 1, 1))
  expect_equal(res$neg_prob, c(1, 1, 0))
  expect_equal(res$lfsr, c(0, 1, 0))
  expect_error(fashr:::compute_posterior_sign_prob(c(1, 2), c(1, -1)),
               "non-negative")
})

test_that("compute_lfsr_summary computes the lfsr of the posterior mixture", {
  fit <- get_test_fit()
  res <- compute_lfsr_summary(fit, index = 2)
  expect_equal(nrow(res), 8)
  expect_true(all(res$lfsr >= 0 & res$lfsr <= 1))
  # Regression test: the lfsr must be the minimum of the mixture-averaged
  # sign probabilities, not the mixture average of per-component minima.
  expect_equal(res$lfsr, pmin(res$pos_prob, res$neg_prob))
  expect_error(compute_lfsr_summary(fit, index = 10), "out of range")
})

test_that("min_lfsr_summary returns one row per dataset with cumulative FSR", {
  fit <- get_test_fit()
  res <- min_lfsr_summary(fit)
  expect_equal(nrow(res), 3)
  expect_equal(res$index, 1:3)
  expect_true(all(res$min_lfsr >= 0 & res$min_lfsr <= 1))
  # fsr is the running mean of the sorted min_lfsr, mapped back to units
  sorted <- res[order(res$min_lfsr), ]
  expect_equal(sorted$fsr, cumsum(sorted$min_lfsr) / seq_len(3))
})
