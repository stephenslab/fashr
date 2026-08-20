test_that("collapse_L uses supplied weights instead of re-estimating", {
  L <- matrix(c(0.5, 0.2, 0.3,
                0.1, 0.6, 0.3), nrow = 2, byrow = TRUE)
  w <- c(0.9, 0.1)
  res <- fashr:::collapse_L(L, log = FALSE, weights = w)
  expect_equal(res$pi_hat_star, w)
  expect_equal(res$L_c[, 1], L[, 1])
  expect_equal(res$L_c[, 2], as.numeric(L[, 2:3] %*% w))
  # weights are normalized when they do not sum to 1
  res2 <- fashr:::collapse_L(L, log = FALSE, weights = c(9, 1))
  expect_equal(res2$pi_hat_star, w)
  # validation
  expect_error(fashr:::collapse_L(L, weights = c(1, 2, 3)), "one entry per")
  expect_error(fashr:::collapse_L(L, weights = c(-1, 2)), "non-negative")
})

test_that("BF_compute reuses estimated prior weights and matches manual computation", {
  grid <- c(0, 0.5, 1)
  L_matrix <- rbind(c(-10, -2, -3),
                    c(-2, -10, -12))
  prior_weight <- c(0.5, 0.3, 0.2)
  obj <- make_mock_fash(L_matrix, grid, prior_weight)
  BF <- BF_compute(obj)
  # Manual computation: row-max rescaled likelihoods, alt weights (0.6, 0.4)
  w_alt <- c(0.3, 0.2) / 0.5
  L_resc <- exp(L_matrix - apply(L_matrix, 1, max))
  BF_manual <- as.numeric(L_resc[, 2:3] %*% w_alt) / L_resc[, 1]
  expect_equal(BF, BF_manual)
})

test_that("BF_compute is robust to strongly negative log-likelihoods (no NaN)", {
  grid <- c(0, 0.5, 1)
  L_matrix <- rbind(c(-5000, -4990, -4995),
                    c(-4990, -5000, -5005))
  obj <- make_mock_fash(L_matrix, grid, c(0.5, 0.3, 0.2))
  BF <- BF_compute(obj)
  expect_true(all(is.finite(BF)))
  expect_gt(BF[1], 1)
  expect_lt(BF[2], 1)
})

test_that("BF_control selects the threshold with the strict 1 + epsilon criterion", {
  # sorted BF = (1, 3): cumulative means are exactly (1, 2).
  # With epsilon = 0 the first crossing is at rank 1 (mu = 1), giving pi0 = 0.5;
  # the strict criterion mu >= 1 + epsilon skips the exact-equality point.
  res_strict <- BF_control(c(3, 1))
  expect_equal(res_strict$pi0_hat_star, 1)
  res_loose <- BF_control(c(3, 1), epsilon = 0)
  expect_equal(res_loose$pi0_hat_star, 0.5)
  # Away from exact equality the default epsilon does not change the result
  set.seed(1)
  BF <- runif(100, 0.2, 5)
  expect_equal(BF_control(BF)$pi0_hat_star, BF_control(BF, epsilon = 0)$pi0_hat_star)
  # pi0 = 1 when the cumulative mean never reaches the threshold
  expect_equal(BF_control(c(0.2, 0.5))$pi0_hat_star, 1)
  expect_error(BF_control(c(1, 2), epsilon = -1), "non-negative")
})

test_that("BF_compute is invariant to the ordering of the PSD grid", {
  grid <- c(0, 0.5, 1)
  L_matrix <- rbind(c(-10, -2, -3), c(-2, -10, -12))
  prior_weight <- c(0.5, 0.3, 0.2)
  obj <- make_mock_fash(L_matrix, grid, prior_weight)
  # Same object with permuted columns (null component not first)
  perm <- c(2, 1, 3)
  obj_perm <- make_mock_fash(L_matrix[, perm, drop = FALSE], grid[perm],
                             prior_weight[perm])
  expect_equal(BF_compute(obj), BF_compute(obj_perm))
})

test_that("BF_update updates weights using the estimated prior and stores BF/lfdr", {
  set.seed(1)
  grid <- c(0, 0.5, 1)
  # 10 null-favoring and 5 alternative-favoring datasets
  L_matrix <- rbind(
    matrix(rep(c(0, -8, -10), each = 10), nrow = 10) + rnorm(30, sd = 0.1),
    matrix(rep(c(-10, -1, 0), each = 5), nrow = 5) + rnorm(15, sd = 0.1)
  )
  obj <- make_mock_fash(L_matrix, grid, c(0.4, 0.35, 0.25))
  updated <- BF_update(obj, plot = FALSE)
  expect_s3_class(updated, "fash")
  expect_length(updated$BF, nrow(L_matrix))
  expect_true(all(is.finite(updated$BF)))
  expect_equal(sum(updated$prior_weights$prior_weight), 1, tolerance = 1e-8)
  # lfdr is the posterior weight of the null component
  null_idx <- which(updated$prior_weights$psd == 0)
  expect_equal(updated$lfdr, updated$posterior_weights[, null_idx])
  # alternative weights keep the ratio from the stored prior (0.35 : 0.25),
  # rather than being re-estimated by mix-SQP
  alt <- updated$prior_weights[updated$prior_weights$psd > 0, ]
  expect_equal(alt$prior_weight[alt$psd == 0.5] / alt$prior_weight[alt$psd == 1],
               0.35 / 0.25, tolerance = 1e-8)
})

test_that("BF machinery errors cleanly without a null component or with one column", {
  grid <- c(0.5, 1)
  L_matrix <- rbind(c(-1, -2), c(-2, -1))
  obj <- make_mock_fash(L_matrix, grid, c(0.6, 0.4))
  expect_error(BF_compute(obj), "null component")
  obj1 <- make_mock_fash(matrix(c(-1, -2), ncol = 1), 0, 1)
  expect_error(BF_compute(obj1), "at least two columns")
})
