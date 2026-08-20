# Shared small gaussian fit, computed lazily and cached across test files.
test_fit_cache <- new.env(parent = emptyenv())

get_test_fit <- function() {
  if (is.null(test_fit_cache$fit)) {
    set.seed(1)
    datasets <- lapply(1:3, function(i) {
      x <- seq(0, 4, length.out = 8)
      f <- if (i %% 2 == 0) function(x) 1.5 * sin(x) else function(x) 0.2 * x
      data.frame(y = f(x) + rnorm(8, sd = 0.3), x = x, sd = 0.3)
    })
    test_fit_cache$fit <- suppressWarnings(
      fash(Y = "y", smooth_var = "x", S = "sd", data_list = datasets,
           grid = c(0, 0.5, 1), order = 2, num_basis = 10,
           likelihood = "gaussian", verbose = FALSE))
  }
  test_fit_cache$fit
}

# A synthetic fash object (no TMB fit required) for testing BF machinery.
make_mock_fash <- function(L_matrix, grid, prior_weight) {
  n <- nrow(L_matrix)
  posterior_weights <- matrix(1 / length(grid), nrow = n, ncol = length(grid))
  rownames(posterior_weights) <- paste0("Dataset_", seq_len(n))
  colnames(posterior_weights) <- as.character(grid)
  structure(list(
    L_matrix = L_matrix,
    psd_grid = grid,
    prior_weights = data.frame(psd = grid, prior_weight = prior_weight),
    posterior_weights = posterior_weights,
    lfdr = posterior_weights[, 1]
  ), class = "fash")
}
