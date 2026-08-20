#' Collapse a Likelihood Matrix for Bayes Factor Computation
#'
#' This function collapses a likelihood matrix into a 2-column matrix, reweighting the likelihood
#' under the alternative hypothesis. By default the weights of the alternative components are
#' estimated with the mix-SQP algorithm, but pre-estimated weights (e.g., the prior weights
#' already fitted by \code{fash} or \code{BF_update}) can be supplied through \code{weights}.
#'
#' @param L A numeric matrix representing the likelihoods. Rows correspond to datasets, and
#'   columns correspond to mixture components (including the null component in the first column).
#' @param log A logical value. If \code{TRUE}, treats \code{L} as a log-likelihood matrix.
#' @param weights An optional numeric vector of non-negative weights for the alternative
#'   components (one per non-null column of \code{L}, in the same order). When provided, these
#'   weights are normalized and used directly instead of re-estimating them with mix-SQP.
#'   This is preferable when the mixture weights have already been estimated, for example
#'   with a Dirichlet penalty (\code{penalty > 1} in \code{fash}), since re-estimating them
#'   here would silently discard the penalty.
#'
#' @return A list containing:
#' \describe{
#'   \item{L_c}{A 2-column matrix where the first column corresponds to the null likelihood
#'   and the second column corresponds to the reweighted alternative likelihood.}
#'   \item{pi_hat_star}{A numeric vector of mixture weights used for the alternative hypothesis.}
#' }
#'
#' @examples
#' # Example likelihood matrix (log-space)
#' set.seed(1)
#' L <- matrix(abs(rnorm(20)), nrow = 5, ncol = 4)
#' collapse_result <- fashr:::collapse_L(L, log = FALSE)
#' print(collapse_result$L_c)
#'
#' # Using pre-estimated weights for the alternative components
#' collapse_result2 <- fashr:::collapse_L(L, log = FALSE, weights = c(0.2, 0.3, 0.5))
#' print(collapse_result2$L_c)
#'
#' @importFrom mixsqp mixsqp
#'
#' @keywords internal
#'
collapse_L <- function(L, log = FALSE, weights = NULL) {
  if (!is.null(weights)) {
    if (length(weights) != ncol(L) - 1) {
      stop("weights must have one entry per non-null (alternative) column of L.")
    }
    if (any(weights < 0)) {
      stop("weights must be non-negative.")
    }
    if (sum(weights) > 0) {
      pi_hat_star <- weights / sum(weights)
    } else {
      # No prior mass on the alternative components; keep the zero weights so
      # the collapsed alternative likelihood is identically zero.
      pi_hat_star <- weights
    }
    if (log) {
      L <- exp(L - apply(L, 1, max))
    }
  } else if (ncol(L) > 1) {
    pi_hat_star <- mixsqp::mixsqp(L = L,
                                  log = log,
                                  control = list(verbose = FALSE))$x[-1]
    pi_hat_star <- pi_hat_star / sum(pi_hat_star)
    if (log) {
      L <- exp(L - apply(L, 1, max))
    }
  } else {
    pi_hat_star <- rep(1, nrow(L))
  }

  L_c <- matrix(0, nrow = nrow(L), ncol = 2)
  L_c[, 1] <- L[, 1]
  L_c[, 2] <- (L[, -1, drop = FALSE] %*% pi_hat_star)

  return(list(L_c = L_c, pi_hat_star = pi_hat_star))
}

#' Extract the Collapsed Likelihood Matrix from a FASH Object
#'
#' Internal helper shared by \code{BF_compute} and \code{BF_update}. It rescales each row
#' of the log-likelihood matrix by its maximum before exponentiating (to avoid numerical
#' underflow; Bayes factors and mixture weight estimates are invariant to row scaling),
#' moves the null (PSD = 0) column first, and collapses the alternative columns using the
#' prior weights already estimated in the \code{fash} object.
#'
#' @param fash A \code{fash} object containing the fitted model and likelihood matrix.
#'
#' @return A list with components \code{L_c} (the collapsed 2-column likelihood matrix),
#'   \code{pi_hat_star} (the normalized weights of the alternative components), and
#'   \code{grid_ordered} (the PSD grid with the null component first).
#'
#' @keywords internal
#'
collapse_L_fash <- function(fash) {
  if (ncol(fash$L_matrix) < 2) {
    stop("The likelihood matrix should have at least two columns (one for the null and one for the alternative). Please check your model specification.")
  }

  grid <- fash$psd_grid
  null_col <- which(grid == 0)
  if (length(null_col) != 1) {
    stop("The PSD grid must contain exactly one null component (PSD = 0) to compute Bayes Factors.")
  }

  # Reorder columns so that the null component comes first
  ord <- c(null_col, setdiff(seq_along(grid), null_col))
  grid_ordered <- grid[ord]

  # Rescale each row by its maximum before exponentiating to avoid underflow;
  # Bayes factors are ratios within each row and are invariant to row scaling.
  L <- exp(fash$L_matrix - apply(fash$L_matrix, 1, max))
  L <- L[, ord, drop = FALSE]

  # Reuse the prior weights already estimated in the fash object (possibly
  # with a Dirichlet penalty) rather than re-estimating them here.
  w_full <- numeric(length(grid_ordered))
  w_full[match(fash$prior_weights$psd, grid_ordered)] <- fash$prior_weights$prior_weight
  collapse_result <- collapse_L(L, log = FALSE, weights = w_full[-1])

  return(list(L_c = collapse_result$L_c,
              pi_hat_star = collapse_result$pi_hat_star,
              grid_ordered = grid_ordered))
}

#' Compute Bayes Factors for Each Dataset in a FASH Object
#'
#' This function computes Bayes Factors (BF) for each dataset in a \code{fash} object.
#' The BF is calculated as the ratio of likelihood under the alternative hypothesis
#' to the likelihood under the null hypothesis, where the alternative components are
#' weighted by the prior weights already estimated in the \code{fash} object.
#'
#' @param fash A \code{fash} object containing the fitted model and likelihood matrix.
#'
#' @return A numeric vector of Bayes Factors, where each entry corresponds to a dataset.
#'   A value of \code{Inf} indicates that the null likelihood underflowed relative to the
#'   best-fitting component (overwhelming evidence against the null).
#'
#' @examples
#' set.seed(1)
#' data_list <- list(
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0),
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0)
#' )
#' grid <- seq(0, 2, length.out = 10)
#' fash_obj <- fash(data_list = data_list, Y = "y", smooth_var = "x", grid = grid, likelihood = "poisson", verbose = TRUE)
#'
#' # Compute Bayes Factors
#' BF_values <- BF_compute(fash_obj)
#' print(BF_values)
#'
#' @export
#'
BF_compute <- function(fash){
  L_c <- collapse_L_fash(fash)$L_c
  BF <- L_c[, 2] / L_c[, 1]
  return(BF)
}

#' Perform Bayes Factor-Based Control for Estimating \eqn{\pi_0}
#'
#' This function estimates \eqn{\pi_0}, the proportion of datasets that follow the null hypothesis,
#' using Bayes Factor (BF) control.
#'
#' @param BF A numeric vector of Bayes Factors computed from `BF_compute()`.
#' @param plot A logical value. If \code{TRUE}, generates diagnostic plots for BF control.
#' @param epsilon A small non-negative numeric value making the threshold selection strict:
#'   the BF threshold \eqn{c^*} is the smallest \eqn{c} such that
#'   \eqn{E(BF \mid BF \le c) \ge 1 + \epsilon}. Defaults to
#'   \code{.Machine$double.eps}, which leaves the results essentially unchanged
#'   while ruling out selection at exact equality.
#'
#' @return A list containing:
#' \describe{
#'   \item{mu}{Cumulative mean of sorted Bayes Factors.}
#'   \item{pi0_hat}{Estimated \eqn{\pi_0} values for each BF threshold.}
#'   \item{pi0_hat_star}{Final estimated \eqn{\pi_0} based on the first BF threshold where \eqn{E(BF \mid BF \le c) \ge 1 + \epsilon}.}
#' }
#'
#' @examples
#' set.seed(1)
#' BF_values <- runif(100, 0.5, 5)  # Example Bayes Factors
#' BF_control_results <- BF_control(BF_values, plot = TRUE)
#' print(BF_control_results$pi0_hat_star)
#'
#' @importFrom graphics par
#' @importFrom graphics hist
#' @importFrom graphics abline
#'
#' @export
#'
BF_control <- function(BF, plot = FALSE, epsilon = .Machine$double.eps) {

  if (!is.numeric(epsilon) || length(epsilon) != 1 || epsilon < 0) {
    stop("epsilon must be a single non-negative numeric value.")
  }

  # check if BF is all NA or NaN
  if (all(is.na(BF)) || all(is.nan(BF))) {
    stop("Bayes Factors contain only NA or NaN values. Please consider refitting the model or checking the data if you wish to use the BF-based correction for prior.")
  }

  # if BF contains NA or NaN, provide a warning
  if (any(is.na(BF)) || any(is.nan(BF))) {
    warning("Bayes Factors contain NA or NaN values. These will be ignored in the analysis.")
    BF <- BF[!is.na(BF) & !is.nan(BF)]
  }

  BF_sorted <- sort(BF, decreasing = FALSE)

  mu <- cumsum(BF_sorted) / seq_along(BF_sorted)
  pi0_hat <- seq_along(BF_sorted) / length(BF_sorted)

  # Select the smallest threshold c such that E(BF | BF <= c) >= 1 + epsilon
  crossing <- which(mu >= 1 + epsilon)[1]
  pi0_hat_star <- if (is.na(crossing)) 1 else pi0_hat[crossing]

  if (plot) {
    par(mfrow = c(1, 2))
    hist(log(BF_sorted[is.finite(BF_sorted)]), breaks = 100, freq = TRUE,
         xlab = "log-BF", main = "Histogram of log-BF")  # Avoid log(Inf) in plot
    if (!is.na(crossing)) {
      abline(v = log(BF_sorted[crossing]), col = "red")
    }

    plot(pi0_hat, mu, type = "l", xlab = "est pi0", ylab = "E(BF | BF <= c)", xlim = c(0,1), ylim = c(0,3))
    abline(h = 1 + epsilon, col = "red")
    par(mfrow = c(1, 1))
  }

  return(list(mu = mu, pi0_hat = pi0_hat, pi0_hat_star = pi0_hat_star))
}


#' Update Prior and Posterior Weights Given \eqn{\pi_0} and \eqn{\pi_{alt}}
#'
#' This function updates the prior and posterior weights in a FASH model using the estimated
#' proportion of null datasets (\eqn{\pi_0}) and the reweighted prior under the alternative hypothesis (\eqn{\pi_{alt}}).
#'
#' @param L_matrix A numeric matrix representing the log-likelihoods of datasets across mixture components.
#'   Rows correspond to datasets, and columns correspond to mixture components.
#' @param pi0 A numeric scalar representing the estimated proportion of null datasets.
#'
#' @param pi_alt A numeric vector representing the estimated weights of the alternative components.
#'
#' @param grid A numeric vector representing the grid of Predictive Standard Deviation (PSD) values.
#'
#' @param null_col An integer giving the column of \code{L_matrix} (and position in \code{grid})
#'   corresponding to the null component (PSD = 0). Defaults to 1.
#'
#' @return A list containing:
#' \describe{
#'   \item{prior_weight}{A data frame with two columns:
#'   \describe{
#'      \item{psd}{A numeric vector of PSD values corresponding to non-trivial weights.}
#'      \item{prior_weight}{A numeric vector of prior weights corresponding to the PSD values.}
#'   }}
#'   \item{posterior_weight}{A numeric matrix of posterior weights, where rows correspond to datasets
#'     and columns correspond to non-trivial mixture components.}
#' }
#'
#' @examples
#' # Example usage:
#' set.seed(1)
#' L_matrix <- matrix(rnorm(50), nrow = 10, ncol = 5)
#' pi0_hat <- 0.8
#' pi_alt <- rep(0.25, 4)  # Alternative weights
#' grid <- seq(0, 2, length.out = 5)
#' update_result <- fashr:::fash_prior_posterior_update(L_matrix, pi0_hat, pi_alt, grid)
#'
#' # View updated prior weights
#' print(update_result$prior_weight)
#'
#' # View updated posterior weights
#' print(update_result$posterior_weight)
#'
#' @keywords internal
#'
fash_prior_posterior_update <- function (L_matrix, pi0, pi_alt, grid, null_col = 1) {
  num_datasets <- nrow(L_matrix)
  num_components <- ncol(L_matrix)

  result_weight <- numeric(num_components)
  result_weight[null_col] <- pi0
  result_weight[-null_col] <- pi_alt * (1 - pi0)
  non_trivial <- which(result_weight > 0)

  prior_weight <- data.frame(
    psd = grid[non_trivial],
    prior_weight = result_weight[non_trivial]
  )

  # Compute posterior weights for each dataset
  posterior_weight <- matrix(0, nrow = num_datasets, ncol = length(non_trivial))
  for (i in 1:num_datasets) {
    exp_values <- exp(L_matrix[i, ] - max(L_matrix[i, ]) + log(result_weight))
    normalized_values <- exp_values[non_trivial] / sum(exp_values[non_trivial])
    posterior_weight[i, ] <- normalized_values
  }
  colnames(posterior_weight) <- as.character(grid[non_trivial])
  rownames(posterior_weight) <- rownames(L_matrix)
  # Return results
  return(list(
    prior_weight = prior_weight,
    posterior_weight = posterior_weight
  ))
}

#' Update Prior and Posterior Weights in a FASH Object Using Bayes Factor Control
#'
#' This function updates the prior and posterior weights in a fitted \code{fash} object using
#' Bayes Factor (BF) control. It automatically computes the Bayes Factor (BF), estimates
#' the proportion of null datasets (\eqn{\pi_0}), and updates the model accordingly.
#'
#' @param fash A \code{fash} object containing the fitted model and likelihood matrix.
#' @param plot A logical value. If \code{TRUE}, generates diagnostic plots for BF control.
#' @param epsilon A small non-negative numeric value passed to \code{BF_control},
#'   making the threshold selection strict: the BF threshold \eqn{c^*} is the
#'   smallest \eqn{c} such that \eqn{E(BF \mid BF \le c) \ge 1 + \epsilon}.
#'   Defaults to \code{.Machine$double.eps}, which leaves the results
#'   essentially unchanged while ruling out selection at exact equality.
#'
#' @return The updated \code{fash} object with the following components updated:
#' \describe{
#'   \item{prior_weights}{Updated prior mixture weights reflecting the estimated \eqn{\pi_0}.}
#'   \item{posterior_weights}{Updated posterior mixture weights for each dataset.}
#'   \item{BF}{Computed Bayes Factors for each dataset.}
#'   \item{lfdr}{Local False Discovery Rate (LFDR), the posterior weight of the null component.}
#' }
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item \bold{Computes Bayes Factors (BF)}: The BF is calculated as the ratio of likelihood under
#'         the alternative hypothesis to the likelihood under the null hypothesis. The weights of
#'         the alternative components are taken from the prior weights already estimated in the
#'         \code{fash} object (so that, e.g., a Dirichlet penalty specified via \code{penalty > 1}
#'         is respected), rather than being re-estimated.
#'   \item \bold{Estimates \eqn{\pi_0}}: The function applies BF-based control to estimate
#'         the proportion of null datasets.
#'   \item \bold{Updates prior weights}: The function updates the prior mixture weights to reflect
#'         the estimated null proportion.
#'   \item \bold{Updates posterior weights}: The posterior weights are updated based on the
#'         reweighted prior and likelihood matrix.
#'   \item \bold{Stores the computed Bayes Factors and LFDR}: The function now saves the computed BF
#'         and LFDR in the \code{fash} object for further analysis.
#' }
#'
#' @examples
#'
#' # Example usage:
#' set.seed(1)
#' data_list <- list(
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0),
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0)
#' )
#' grid <- seq(0, 2, length.out = 10)
#' fash_obj <- fash(data_list = data_list, Y = "y", smooth_var = "x",
#'                  grid = grid, likelihood = "poisson", verbose = TRUE)
#'
#' # Update prior and posterior weights using BF control
#' fash_updated <- BF_update(fash_obj, plot = TRUE)
#'
#' # Access updated components
#' print(fash_updated$prior_weights)
#' print(fash_updated$posterior_weights)
#' print(fash_updated$BF)
#' print(fash_updated$lfdr)
#'
#' @export
#'
BF_update <- function (fash, plot = FALSE, epsilon = .Machine$double.eps) {

  # Collapse the likelihood matrix once, reusing the estimated prior weights
  collapse_result <- collapse_L_fash(fash)
  L_c <- collapse_result$L_c
  pi_alt <- collapse_result$pi_hat_star
  grid_ordered <- collapse_result$grid_ordered

  # Compute Bayes Factors
  BF <- L_c[, 2] / L_c[, 1]

  # check if BF is all NA or NaN
  if (all(is.na(BF)) || all(is.nan(BF))) {
    # provide a warning and return the fash object without updating
    warning("Bayes Factors contain only NA or NaN values. BF-based correction cannot be applied. Returning the original fash object without updates. Please consider refitting the model or checking the data if you wish to use the BF-based correction for prior.")
    return(fash)
  }

  # Perform BF control
  BF_res <- BF_control(BF, plot = plot, epsilon = epsilon)
  pi0_hat <- BF_res$pi0_hat_star

  # Reorder the likelihood matrix to match grid_ordered (null component first)
  null_col <- which(fash$psd_grid == 0)
  ord <- c(null_col, setdiff(seq_along(fash$psd_grid), null_col))
  L_matrix <- fash$L_matrix[, ord, drop = FALSE]
  rownames(L_matrix) <- rownames(fash$posterior_weights)

  # Update prior and posterior weights
  update_res <- fash_prior_posterior_update(L_matrix = L_matrix,
                  pi0 = pi0_hat, pi_alt = pi_alt, grid = grid_ordered,
                  null_col = 1)

  # Update fash object
  fash$prior_weights <- update_res$prior_weight
  fash$posterior_weights <- update_res$posterior_weight
  fash$BF <- BF
  null_idx <- which(update_res$prior_weight$psd == 0)
  if (length(null_idx) == 1) {
    fash$lfdr <- fash$posterior_weights[, null_idx]
  } else {
    warning("The updated prior weight of the null component (PSD = 0) is zero; lfdr is set to 0 for all datasets.")
    fash$lfdr <- rep(0, nrow(fash$posterior_weights))
  }

  return(fash)
}
