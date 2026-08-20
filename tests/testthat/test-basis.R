test_that("global_poly_helper builds the polynomial design matrix", {
  x <- c(0, 0.5, 1, 2)
  X <- fashr:::global_poly_helper(x, p = 3)
  expect_equal(dim(X), c(4, 3))
  expect_equal(X[, 1], rep(1, 4))
  expect_equal(X[, 2], x)
  expect_equal(X[, 3], x^2)
})

test_that("local_poly_helper builds the O-spline design matrix", {
  knots <- c(0, 1, 2, 3)
  x <- seq(0, 3, by = 0.5)
  p <- 2
  B <- fashr:::local_poly_helper(knots = knots, refined_x = x, p = p)
  expect_equal(dim(B), c(length(x), length(knots) - 1))
  # Basis functions vanish at and before their starting knot
  expect_equal(B[x <= 1, 2], rep(0, sum(x <= 1)))
  # Inside the knot interval the basis is (x - knot)^p / p!
  expect_equal(B[x == 0.5, 1], 0.5^p / factorial(p))
  # Beyond the interval the basis continues as a polynomial of degree p - 1
  # with matching derivatives at the right knot: value + slope extension
  expect_equal(B[x == 2, 1], 1 / factorial(2) + (2 - 1) * 1)
})

test_that("compute_weights_precision_helper returns diag of knot spacings", {
  x <- c(0, 0.2, 0.5, 1)
  P <- fashr:::compute_weights_precision_helper(x)
  expect_equal(P, diag(diff(x)))
})
