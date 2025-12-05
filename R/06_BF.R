#' Collapse a Likelihood Matrix for Bayes Factor Computation
#'
#' This function collapses a likelihood matrix into a 2-column matrix, reweighting the likelihood
#' under the alternative hypothesis. This is done using the mix-SQP algorithm to estimate
#' the optimal mixture weights under the alternative hypothesis.
#'
#' @param L A numeric matrix representing the likelihoods. Rows correspond to datasets, and
#'   columns correspond to mixture components (including the null component in the first column).
#' @param log A logical value. If \code{TRUE}, treats \code{L} as a log-likelihood matrix.
#'
#' @return A list containing:
#' \describe{
#'   \item{L_c}{A 2-column matrix where the first column corresponds to the null likelihood
#'   and the second column corresponds to the reweighted alternative likelihood.}
#'   \item{pi_hat_star}{A numeric vector of mixture weights estimated under the alternative hypothesis.}
#' }
#'
#' @examples
#' # Example likelihood matrix (log-space)
#' set.seed(1)
#' L <- matrix(abs(rnorm(20)), nrow = 5, ncol = 4)
#' collapse_result <- fashr:::collapse_L(L, log = FALSE)
#' print(collapse_result$L_c)
#'
#' @importFrom mixsqp mixsqp
#'
#' @keywords internal
#'
collapse_L <- function(L, log = FALSE) {
  if (ncol(L) > 1) {
    pi_hat_star <- mixsqp::mixsqp(L = L,
                                  log = log,
                                  control = list(verbose = FALSE))$x[-1]
    pi_hat_star <- pi_hat_star / sum(pi_hat_star)
  } else {
    pi_hat_star <- rep(1, nrow(L))
  }

  L_c <- matrix(0, nrow = nrow(L), ncol = 2)
  L_c[, 1] <- L[, 1]
  L_c[, 2] <- (L[, -1, drop = FALSE] %*% pi_hat_star)

  return(list(L_c = L_c, pi_hat_star = pi_hat_star))
}

#' Compute Bayes Factors for Each Dataset in a FASH Object
#'
#' This function computes Bayes Factors (BF) for each dataset in a \code{fash} object.
#' The BF is calculated as the ratio of likelihood under the alternative hypothesis
#' to the likelihood under the null hypothesis.
#'
#' @param fash A \code{fash} object containing the fitted model and likelihood matrix.
#'
#' @return A numeric vector of Bayes Factors, where each entry corresponds to a dataset.
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
  # Check the number of columns in L_matrix
  if (ncol(fash$L_matrix) < 2) {
    stop("The likelihood matrix should have at least two columns (one for the null and one for the alternative). Please check your model specification.")
  }

  L <- exp(fash$L_matrix)
  L_c <- collapse_L(L, log = FALSE)$L_c
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
#'
#' @return A list containing:
#' \describe{
#'   \item{mu}{Cumulative mean of sorted Bayes Factors.}
#'   \item{pi0_hat}{Estimated \eqn{\pi_0} values for each BF threshold.}
#'   \item{pi0_hat_star}{Final estimated \eqn{\pi_0} based on the first BF threshold where \eqn{E(BF) \geq 1}.}
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
BF_control <- function(BF, plot = FALSE) {

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

  pi0_hat_star <- if (max(mu, na.rm = TRUE) < 1) 1 else pi0_hat[which(mu >= 1)[1]]

  if (plot) {
    par(mfrow = c(1, 2))
    hist(log(BF_sorted[is.finite(BF_sorted)]), breaks = 100, freq = TRUE,
         xlab = "log-BF", main = "Histogram of log-BF")  # Avoid log(Inf) in plot
    abline(v = log(BF_sorted[which(mu >= 1)[1]]), col = "red")

    plot(pi0_hat, mu, type = "l", xlab = "est pi0", ylab = "E(BF | BF <= c)", xlim = c(0,1), ylim = c(0,3))
    abline(h = 1, col = "red")
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
#' pi_alt <- rep(0.2, 4)  # Alternative weights
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
fash_prior_posterior_update <- function (L_matrix, pi0, pi_alt, grid) {
  num_datasets <- nrow(L_matrix)
  num_components <- ncol(L_matrix)

  result_weight <- c(pi0, pi_alt * (1 - pi0))
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
#'
#' @return The updated \code{fash} object with the following components updated:
#' \describe{
#'   \item{prior_weights}{Updated prior mixture weights reflecting the estimated \eqn{\pi_0}.}
#'   \item{posterior_weights}{Updated posterior mixture weights for each dataset.}
#'   \item{BF}{Computed Bayes Factors for each dataset.}
#'   \item{lfdr}{Local False Discovery Rate (LFDR), extracted as the first column of `posterior_weights`.}
#' }
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item \bold{Computes Bayes Factors (BF)}: The BF is calculated as the ratio of likelihood under
#'         the alternative hypothesis to the likelihood under the null hypothesis.
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
BF_update <- function (fash, plot = FALSE) {

  # Check the number of columns in L_matrix
  if (ncol(fash$L_matrix) < 2) {
    stop("The likelihood matrix should have at least two columns (one for the null and one for the alternative). Please check your model specification.")
  }

  # Compute Lc
  L <- exp(fash$L_matrix)
  L_c <- collapse_L(L, log = FALSE)$L_c
  pi_alt <- collapse_L(L, log = FALSE)$pi_hat_star

  # Compute Bayes Factors
  BF <- L_c[, 2] / L_c[, 1]

  # check if BF is all NA or NaN
  if (all(is.na(BF)) || all(is.nan(BF))) {
    # provide a warning and return the fash object without updating
    warning("Bayes Factors contain only NA or NaN values. BF-based correction cannot be applied. Returning the original fash object without updates. Please consider refitting the model or checking the data if you wish to use the BF-based correction for prior.")
    return(fash)
  }

  # Perform BF control
  BF_res <- BF_control(BF, plot = plot)
  pi0_hat <- BF_res$pi0_hat_star

  L_matrix <- fash$L_matrix
  rownames(L_matrix) <- rownames(fash$posterior_weights)

  # Update prior and posterior weights
  update_res <- fash_prior_posterior_update(L_matrix = L_matrix,
                  pi0 = pi0_hat, pi_alt = pi_alt, grid = fash$psd_grid)

  # Update fash object
  fash$prior_weights <- update_res$prior_weight
  fash$posterior_weights <- update_res$posterior_weight
  fash$BF <- BF
  fash$lfdr <- fash$posterior_weights[, 1]  # LFDR is the first column of posterior_weights

  return(fash)
}
