test_that("predict.fash returns summaries and samples of the right shape", {
  fit <- get_test_fit()
  set.seed(1)
  summ <- predict(fit, index = 1, M = 50)
  expect_s3_class(summ, "data.frame")
  expect_equal(nrow(summ), 8)
  expect_named(summ, c("x", "mean", "median", "lower", "upper"))
  expect_true(all(summ$lower <= summ$upper))
  samps <- predict(fit, index = 1, smooth_var = seq(0, 4, length.out = 11),
                   only.samples = TRUE, M = 30)
  expect_equal(dim(samps), c(11, 30))
  expect_true(all(is.finite(samps)))
  expect_error(predict(fit, index = 4), "out of range")
})

test_that("predict.fash works when the prior collapses to a single PSD component", {
  fit <- get_test_fit()
  # Force a degenerate prior with a single non-trivial PSD value larger than 1;
  # regression test for the sample() length-1 shorthand.
  fit$prior_weights <- data.frame(psd = 2, prior_weight = 1)
  fit$posterior_weights <- matrix(1, nrow = 3, ncol = 1,
                                  dimnames = list(rownames(fit$posterior_weights), "2"))
  set.seed(1)
  samps <- predict(fit, index = 1, only.samples = TRUE, M = 10)
  expect_equal(dim(samps), c(8, 10))
  expect_true(all(is.finite(samps)))
})

test_that("predict.fash works when the prior collapses to the null component", {
  fit <- get_test_fit()
  fit$prior_weights <- data.frame(psd = 0, prior_weight = 1)
  fit$posterior_weights <- matrix(1, nrow = 3, ncol = 1,
                                  dimnames = list(rownames(fit$posterior_weights), "0"))
  set.seed(1)
  samps <- predict(fit, index = 1, only.samples = TRUE, M = 10)
  expect_equal(dim(samps), c(8, 10))
  expect_true(all(is.finite(samps)))
})

test_that("compute_lfsr_sampling returns probabilities of the right shape", {
  fit <- get_test_fit()
  set.seed(1)
  res <- compute_lfsr_sampling(fit, index = 2, M = 100)
  expect_length(res$lfsr, 8)
  expect_true(all(res$lfsr >= 0 & res$lfsr <= 1))
  # samples are re-centered at the first point, so its lfsr is degenerate
  expect_equal(res$lfsr[1], 1)
  expect_equal(res$lfsr, pmin(res$pos_prob, res$neg_prob))
})
