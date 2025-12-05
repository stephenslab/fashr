#' Compute Local False Sign Rate (LFSR) for a Single Dataset
#'
#' This function computes the local false sign rate (LFSR) for a specific dataset within a \code{fash} object
#' using posterior sampling. It estimates the probability that the sign of the effect is positive or negative
#' at each grid point and returns the full probability vectors.
#'
#' @param fash_fit A \code{fash} object containing posterior samples.
#' @param index An integer specifying the dataset index for which to compute the LFSR.
#' @param smooth_var A numeric vector specifying refined x values for prediction.
#'   If \code{NULL}, defaults to the dataset's original x values.
#' @param M An integer specifying the number of posterior samples to generate.
#' @param deriv An integer specifying the order of the derivative to compute.
#'
#' @return A list containing:
#' \describe{
#'   \item{lfsr}{A numeric vector of LFSR values, where each entry corresponds to a grid point in \code{smooth_var}.}
#'   \item{pos_prob}{A numeric vector of probabilities that the effect is positive at each grid point.}
#'   \item{neg_prob}{A numeric vector of probabilities that the effect is negative at each grid point.}
#' }
#'
#' @examples
#' # Example fash object (assuming it has been fitted)
#' set.seed(1)
#' data_list <- list(
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0),
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0)
#' )
#' grid <- seq(0, 2, length.out = 10)
#' fash_obj <- fash(data_list = data_list, Y = "y", smooth_var = "x", grid = grid, likelihood = "poisson", verbose = TRUE)
#'
#' # Compute LFSR for a single dataset
#' lfsr_result <- compute_lfsr_sampling(fash_obj, index = 1, smooth_var = seq(0, 5, by = 0.1), M = 3000)
#' print(lfsr_result$lfsr)  # Print the LFSR values
#' print(lfsr_result$pos_prob)  # Print the positive probability
#' print(lfsr_result$neg_prob)  # Print the negative probability
#'
#' @importFrom stats pnorm
#'
#' @export
#'
compute_lfsr_sampling <- function(fash_fit, index, smooth_var = NULL, M = 3000, deriv = 0) {
  # Validate input
  if (!inherits(fash_fit, "fash")) stop("fash_fit must be a `fash` object.")
  if (!is.numeric(index) || index <= 0 || index > length(fash_fit$fash_data$data_list))
    stop("index must be a valid dataset index.")
  if (!is.numeric(M) || M <= 0 || M %% 1 != 0) stop("M must be a positive integer.")

  # Obtain posterior samples for the specified dataset
  sample_i <- predict(fash_fit, index = index, smooth_var = smooth_var, only.samples = TRUE, M = M, deriv = deriv)

  # Compute LFSR
  sample_i_tilde <- apply(sample_i, 2, function(x) x - x[1])  # Re-center samples
  pos_prob <- rowMeans(sample_i_tilde >= 0)  # Faster than apply()
  neg_prob <- rowMeans(sample_i_tilde <= 0)
  lfsr <- pmin(pos_prob, neg_prob)

  return(list(
    lfsr = lfsr,
    pos_prob = pos_prob,
    neg_prob = neg_prob
  ))
}

#' Compute Minimum Local False Sign Rate (LFSR) from Posterior Samples
#'
#' This function computes the minimum local false sign rate (LFSR) for each dataset in a \code{fash} object.
#' It estimates the probability that the sign of the effect is positive or negative at each x value
#' and returns a data frame.
#'
#' @param fash_fit A \code{fash} object containing posterior samples.
#' @param smooth_var A numeric vector specifying refined x values for prediction.
#'   If \code{NULL}, defaults to the dataset's original x values.
#' @param M An integer specifying the number of posterior samples to generate.
#' @param num_cores An integer specifying the number of cores to use for parallel processing.
#' @param deriv An integer specifying the order of the derivative to compute.
#'
#'
#' @return A data frame containing:
#'
#' \describe{
#'   \item{index}{The dataset index.}
#'   \item{min_lfsr}{The minimum LFSR computed for each dataset.}
#'   \item{fsr}{The cumulative false sign rate (FSR).}
#' }
#'
#' @examples
#' # Example fash object (assuming it has been fitted)
#' set.seed(1)
#' data_list <- list(
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0),
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0)
#' )
#' grid <- seq(0, 2, length.out = 10)
#' fash_obj <- fash(data_list = data_list, Y = "y", smooth_var = "x", grid = grid, likelihood = "poisson", verbose = TRUE)
#'
#' # Compute min LFSR with sequential execution
#' result <- min_lfsr_sampling(fash_obj, num_cores = 1)
#'
#' # Compute min LFSR with parallel execution
#' result_parallel <- min_lfsr_sampling(fash_obj, num_cores = 2)
#'
#' @importFrom parallel mclapply
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
min_lfsr_sampling <- function(fash_fit, smooth_var = NULL, M = 3000, num_cores = 1, deriv = 0) {
  datasets <- fash_fit$fash_data$data_list
  n_datasets <- length(datasets)

  # Validate input
  if (!inherits(fash_fit, "fash")) stop("fash_fit must be a `fash` object.")
  if (!is.numeric(M) || M <= 0 || M %% 1 != 0) stop("M must be a positive integer.")
  if (!is.numeric(num_cores) || num_cores < 1 || num_cores %% 1 != 0) stop("num_cores must be a positive integer.")

  # Define the function to compute min LFSR for a single dataset
  compute_min_lfsr <- function(i) {
    sample_i <- predict(fash_fit, index = i, smooth_var = smooth_var, only.samples = TRUE, M = M, deriv = deriv)
    sample_i_tilde <- apply(sample_i, 2, function(x) x - x[1])  # Re-center samples
    pos_prob <- rowMeans(sample_i_tilde >= 0)  # Faster than apply()
    neg_prob <- rowMeans(sample_i_tilde <= 0)
    lfsr <- pmin(pos_prob, neg_prob)
    return(c(i, min(lfsr)))
  }

  # Parallel execution
  if (num_cores > 1) {
    results_list <- parallel::mclapply(1:n_datasets, compute_min_lfsr, mc.cores = num_cores)
  } else {
    # Sequential execution with progress bar
    pb <- utils::txtProgressBar(min = 0, max = n_datasets, style = 3)
    results_list <- lapply(1:n_datasets, function(i) {
      utils::setTxtProgressBar(pb, i)
      compute_min_lfsr(i)
    })
    close(pb)
  }

  # Convert results to a data frame
  results_mat <- do.call(rbind, results_list)
  lfsr_df <- data.frame(index = results_mat[, 1], min_lfsr = results_mat[, 2])
  lfsr_df <- lfsr_df[order(lfsr_df$min_lfsr), ]

  # Compute cumulative false sign rate (FSR)
  lfsr_df$fsr <- cumsum(lfsr_df$min_lfsr) / seq_len(n_datasets)

  # reorder back to original index
  lfsr_df <- lfsr_df[order(lfsr_df$index), ]

  return(lfsr_df)
}

#' Compute Probability of Being Positive or Negative
#'
#' Given posterior means and variances, this function computes the probability
#' that the posterior distribution is positive or negative at each point.
#'
#' @keywords internal
#'
compute_posterior_sign_prob <- function(mu, sigma2) {
  # Validate inputs
  if (length(mu) != length(sigma2)) stop("mu and sigma2 must have the same length.")
  if (any(sigma2 < 0)) stop("sigma2 must be non-negative.")

  # Convert variance to standard deviation
  sigma <- sqrt(sigma2)

  # Initialize probability vectors
  pos_prob <- numeric(length(mu))
  neg_prob <- numeric(length(mu))

  # Case when sigma = 0 (point mass)
  zero_sigma_idx <- sigma == 0
  nonzero_sigma_idx <- !zero_sigma_idx

  # Compute probabilities using normal CDF for nonzero sigma
  if (any(nonzero_sigma_idx)) {
    pos_prob[nonzero_sigma_idx] <- 1 - pnorm(-mu[nonzero_sigma_idx] / sigma[nonzero_sigma_idx])
    neg_prob[nonzero_sigma_idx] <- pnorm(-mu[nonzero_sigma_idx] / sigma[nonzero_sigma_idx])
  }

  # Handle cases where sigma = 0 (point mass at mu)
  if (any(zero_sigma_idx)) {
    pos_prob[zero_sigma_idx] <- ifelse(mu[zero_sigma_idx] >= 0, 1, 0)
    neg_prob[zero_sigma_idx] <- ifelse(mu[zero_sigma_idx] <= 0, 1, 0)
  }

  # Compute Local False Sign Rate (LFSR)
  lfsr <- pmin(pos_prob, neg_prob)

  # Return a structured data frame
  return(data.frame(pos_prob = pos_prob, neg_prob = neg_prob, lfsr = lfsr))
}

#' Compute Marginal Mean and Variance for a Single PSD Value
#'
#' Computes the posterior mean and variance for a given dataset, refined x-values,
#' and a specific PSD value.
#'
#' @importFrom TMB MakeADFun
#' @importFrom numDeriv jacobian
#' @importFrom Matrix forceSymmetric
#' @importFrom stats nlminb
#'
#' @keywords internal
#'
compute_marginal_mean_var_once <- function(data_i, refined_x, psd_iwp, Si = NULL, Omegai = NULL, num_basis = 30, betaprec = 1e-6, order = 2, pred_step = 1, likelihood, deriv = 0) {
  # Create the tmbdat object using the existing helper function
  tmbdat <- fash_set_tmbdat(data_i, Si, Omegai, num_basis = num_basis, betaprec = betaprec, order = order)

  # Extract smoothing variables and response
  y <- data_i$y
  x <- data_i$x
  offset <- data_i$offset

  # Generate spline knots for the smoothing variable
  knots <- seq(min(x), max(x), length.out = num_basis)

  # check if deriv is not strictly smaller than order
  if(deriv >= order){
    stop("deriv must be strictly smaller than order.")
  }

  if (psd_iwp != 0) {
    B_refined <- local_poly_helper(knots = knots, refined_x = refined_x, p = (order-deriv))
    tmbdat$sigmaIWP <- psd_iwp / sqrt((pred_step ^ ((2 * order) - 1)) / (((2 * order) - 1) * (factorial(order - 1) ^ 2)))
  } else{
    return(data.frame(mean = rep(0, length(refined_x)), var = rep(0, length(refined_x))))
  }

  # Determine the DLL for TMB based on likelihood and error structure
  DLL <- switch(likelihood,
                "gaussian" = if (!is.null(Si)) {
                  if (psd_iwp != 0) "Gaussian_ind" else "Gaussian_ind_fixed"
                } else {
                  if (psd_iwp != 0) "Gaussian_dep" else "Gaussian_dep_fixed"
                },
                "poisson" = if (psd_iwp != 0) "Poisson_ind" else "Poisson_ind_fixed",
                stop("Unknown likelihood function. Choose 'gaussian' or 'poisson'."))

  # Define initial parameter values
  tmbparams <- list(
    W = rep(0, if (psd_iwp != 0) ncol(tmbdat$X) + ncol(tmbdat$B) else ncol(tmbdat$X))
  )

  # Create the TMB model
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    DLL = DLL,
    silent = TRUE
  )

  # Define Hessian for optimization
  ff$he <- function(w) numDeriv::jacobian(ff$gr, w)

  # Optimize the model
  opt <- stats::nlminb(
    start = ff$par,
    objective = ff$fn,
    gradient = ff$gr,
    hessian = ff$he,
    control = list(eval.max = 20000, iter.max = 20000)
  )

  # Extract precision matrix from Hessian
  prec_matrix <- Matrix::forceSymmetric(ff$he(opt$par))
  var_matrix <- B_refined %*% solve(prec_matrix)[(1:ncol(B_refined)), (1:ncol(B_refined))] %*% t(B_refined)
  var <- diag(var_matrix)
  mean <- B_refined %*% opt$par[1:ncol(B_refined), drop = FALSE]
  data.frame(mean = mean, var = var)
}


#' Output all the marginal mean and variance for each PSD value
#' @keywords internal
compute_marginal_mean_var <- function(data_i, psd_values, refined_x,
                                      Si = NULL, Omegai = NULL, num_basis = 30, betaprec = 1e-6,
                                      order = 2, pred_step = 1, likelihood, deriv = 0) {

  # Initialize matrices to store results
  mean_matrix <- matrix(0, nrow = length(refined_x), ncol = length(psd_values))
  var_matrix <- matrix(0, nrow = length(refined_x), ncol = length(psd_values))

  # Loop through each PSD value and compute mean and variance
  for (j in seq_along(psd_values)) {
    psd_iwp <- psd_values[j]
    mean_var_df <- compute_marginal_mean_var_once(
      data_i, refined_x, psd_iwp,
      Si = Si, Omegai = Omegai, num_basis = num_basis, betaprec = betaprec,
      order = order, pred_step = pred_step, likelihood = likelihood, deriv = deriv
    )

    mean_matrix[, j] <- mean_var_df$mean
    var_matrix[, j] <- mean_var_df$var
  }

  return(list(mean = mean_matrix, var = var_matrix))

}



#' Compute Local False Sign Rate (LFSR) from Marginal Posterior Mean and Variance
#'
#' This function computes the local false sign rate (LFSR) at each evaluation point
#' for a specific dataset in a \code{fash} object using the posterior mean and variance
#' instead of sampling-based methods.
#'
#' @param object A \code{fash} object containing the fitted results.
#' @param index An integer specifying the dataset index.
#' @param smooth_var A numeric vector specifying refined x values for evaluation.
#'   If \code{NULL}, defaults to the dataset's original x values.
#' @param deriv An integer specifying the order of the derivative to compute.
#'
#' @return A data frame with the following columns:
#'
#' \describe{
#'   \item{x}{The refined x values where posterior probabilities are evaluated.}
#'   \item{pos_prob}{The posterior probability that the function is positive at each x.}
#'   \item{neg_prob}{The posterior probability that the function is negative at each x.}
#'   \item{lfsr}{The local false sign rate (LFSR) at each x, computed as \code{pmin(pos_prob, neg_prob)}.}
#' }
#'
#' @examples
#' # Example usage
#' set.seed(1)
#' data_list <- list(
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0),
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0)
#' )
#' grid <- seq(0, 2, length.out = 10)
#' fash_obj <- fash(data_list = data_list, Y = "y", smooth_var = "x", grid = grid, likelihood = "poisson", verbose = TRUE)
#'
#' # Compute LFSR summary for dataset 1
#' lfsr_summary_df <- compute_lfsr_summary(fash_obj, index = 1)
#' print(lfsr_summary_df)
#'
#' @export
compute_lfsr_summary <- function(object, index = 1, smooth_var = NULL, deriv = 0) {
  # Validate input
  if (!inherits(object, "fash")) {
    stop("Input must be a `fash` object.")
  }
  if (index < 1 || index > length(object$posterior_weights)) {
    stop("Index is out of range for the datasets in the `fash` object.")
  }

  # Extract dataset-specific components
  data_i <- object$fash_data$data_list[[index]]
  Si <- object$fash_data$S[[index]]
  Omegai <- object$fash_data$Omega[[index]]
  psd_values <- object$prior_weights$psd
  posterior_weights <- object$posterior_weights[index, ]

  # Use smooth_var if provided; otherwise, default to the dataset's x values
  refined_x <- if (!is.null(smooth_var)) smooth_var else data_i$x

  # Retrieve settings from fash object
  settings <- object$settings

  # Compute marginal mean and variance for each PSD value
  marginal_stats <- compute_marginal_mean_var(
    data_i = data_i,
    psd_values = psd_values,
    refined_x = refined_x,
    Si = Si,
    Omegai = Omegai,
    num_basis = settings$num_basis,
    betaprec = settings$betaprec,
    order = settings$order,
    pred_step = settings$pred_step,
    likelihood = settings$likelihood,
    deriv = deriv
  )

  # Compute probability of being positive/negative for each PSD value
  prob_list <- lapply(seq_along(psd_values), function(j) {
    compute_posterior_sign_prob(marginal_stats$mean[, j], marginal_stats$var[, j])
  })

  # Convert list of data frames into matrices
  pos_prob_matrix <- do.call(cbind, lapply(prob_list, function(df) df$pos_prob))
  neg_prob_matrix <- do.call(cbind, lapply(prob_list, function(df) df$neg_prob))
  lfsr_matrix <- do.call(cbind, lapply(prob_list, function(df) df$lfsr))

  # Compute final weighted probability across all PSD values
  weighted_pos_prob <- as.numeric(pos_prob_matrix %*% posterior_weights)
  weighted_neg_prob <- as.numeric(neg_prob_matrix %*% posterior_weights)
  weighted_lfsr <- as.numeric(lfsr_matrix %*% posterior_weights)

  # Return final results as a data frame
  return(data.frame(
    x = refined_x,
    pos_prob = weighted_pos_prob,
    neg_prob = weighted_neg_prob,
    lfsr = weighted_lfsr
  ))
}


#' Compute Minimum Local False Sign Rate (LFSR) for All Datasets
#'
#' This function computes the minimum LFSR for each dataset in a \code{fash} object using
#' the posterior mean and variance instead of sampling-based methods.
#'
#' @param object A \code{fash} object containing the fitted results.
#' @param smooth_var A numeric vector specifying refined x values for evaluation.
#'   If \code{NULL}, defaults to the dataset's original x values.
#' @param num_cores An integer specifying the number of cores to use for parallel processing.
#' @param deriv An integer specifying the order of the derivative to compute.
#'
#' @return A data frame containing:
#'
#' \describe{
#'   \item{index}{The dataset index.}
#'   \item{min_lfsr}{The minimum LFSR computed for each dataset.}
#'   \item{fsr}{The cumulative false sign rate (FSR).}
#' }
#'
#' @examples
#' # Example usage
#' set.seed(1)
#' data_list <- list(
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0),
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0)
#' )
#' grid <- seq(0, 2, length.out = 10)
#' fash_obj <- fash(data_list = data_list, Y = "y", smooth_var = "x", grid = grid, likelihood = "poisson", verbose = TRUE)
#'
#' # Compute min LFSR for all datasets sequentially
#' result <- min_lfsr_summary(fash_obj, num_cores = 1)
#' print(result)
#'
#' # Compute min LFSR for all datasets in parallel
#' result_parallel <- min_lfsr_summary(fash_obj, num_cores = 2)
#' print(result_parallel)
#'
#' @importFrom parallel mclapply
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
min_lfsr_summary <- function(object, smooth_var = NULL, num_cores = 1, deriv = 0) {
  datasets <- object$fash_data$data_list
  n_datasets <- length(datasets)

  # Validate input
  if (!inherits(object, "fash")) stop("Input must be a `fash` object.")
  if (!is.numeric(num_cores) || num_cores < 1 || num_cores %% 1 != 0) stop("num_cores must be a positive integer.")

  # Define function to compute min LFSR for a single dataset
  compute_min_lfsr_single <- function(i) {
    result <- compute_lfsr_summary(object, index = i, smooth_var = smooth_var, deriv = deriv)
    return(data.frame(index = i, min_lfsr = min(result$lfsr)))
  }

  # Parallel execution
  if (num_cores > 1) {
    results_list <- parallel::mclapply(1:n_datasets, compute_min_lfsr_single, mc.cores = num_cores)
  } else {
    # Sequential execution with progress bar
    pb <- utils::txtProgressBar(min = 0, max = n_datasets, style = 3)
    results_list <- lapply(1:n_datasets, function(i) {
      utils::setTxtProgressBar(pb, i)
      compute_min_lfsr_single(i)
    })
    close(pb)
  }

  # Convert list of data frames into a single data frame
  lfsr_df <- do.call(rbind, results_list)

  # Sort by min_lfsr
  lfsr_df <- lfsr_df[order(lfsr_df$min_lfsr), ]

  # Compute cumulative false sign rate (FSR)
  lfsr_df$fsr <- cumsum(lfsr_df$min_lfsr) / seq_len(n_datasets)

  # reorder back to original index
  lfsr_df <- lfsr_df[order(lfsr_df$index), ]

  return(lfsr_df)
}


