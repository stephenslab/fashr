#' Fit Model Once and Extract Posterior Sample Paths
#'
#' This function fits the model for a given dataset and specified PSD value, and generates posterior sample paths at refined grid points.
#'
#' @param data_i A single dataset extracted from the \code{data_list} component of \code{fash_set_data}. Must be a list containing \code{y}, \code{x}, and \code{offset}.
#' @param refined_x A numeric vector of grid points where posterior sample paths are evaluated.
#' @param M An integer specifying the number of posterior samples to generate.
#' @param psd_iwp A numeric value specifying the PSD value for the Integrated Wiener Process.
#' @param Si A numeric vector representing the standard errors for the dataset (if applicable).
#' @param Omegai A numeric precision matrix for the dataset (if applicable).
#' @param num_basis An integer specifying the number of O-Spline basis functions.
#' @param betaprec A numeric value representing the precision of the fixed effects coefficients.
#' @param order An integer specifying the order of the Integrated Wiener Process (IWP) prior.
#' @param pred_step A numeric value specifying the prediction step size.
#' @param likelihood A character string specifying the likelihood function to use. Options are `gaussian` and `poisson`.
#' @param deriv An integer specifying the order of the derivative to compute.
#'
#' @return A numeric matrix of posterior sample paths, where rows correspond to grid points in \code{refined_x} and columns correspond to the generated samples.
#'
#' @examples
#' # Example usage
#' set.seed(1)
#' Y <- matrix(rpois(20, lambda = 5), nrow = 4, ncol = 10, byrow = TRUE)
#' smooth_var <- matrix(runif(20), nrow = 4, ncol = 10, byrow = TRUE)
#' offset <- 1
#' fash_data <- fash_set_data(Y, smooth_var, offset)
#' refined_x <- seq(0, 10, length.out = 50)
#' samples <- fashr:::fash_fit_once(data_i = fash_data$data_list[[1]], refined_x = refined_x, M = 100, num_basis = 60,
#'                          psd_iwp = 0.5, Si = NULL, Omegai = NULL, likelihood = "poisson")
#'
#'
#' @importFrom TMB MakeADFun
#' @importFrom numDeriv jacobian
#' @importFrom LaplacesDemon rmvnp
#' @importFrom stats nlminb
#' @importFrom Matrix forceSymmetric
#'
#' @keywords internal
#'
fash_fit_once <- function(data_i, refined_x, M, psd_iwp, Si = NULL, Omegai = NULL, num_basis = 30, betaprec = 1e-6, order = 1, pred_step = 1, likelihood, deriv = 0) {

  # return error if deriv is not strictly smaller than order
  if(deriv >= order){
    stop("deriv must be strictly smaller than order.")
  }

  # Create the tmbdat object using the existing helper function
  tmbdat <- fash_set_tmbdat(data_i, Si, Omegai, num_basis = num_basis, betaprec = betaprec, order = order)

  # Extract smoothing variables and response
  y <- data_i$y
  x <- data_i$x
  offset <- data_i$offset

  # Generate spline knots for the smoothing variable
  knots <- seq(min(x), max(x), length.out = num_basis)

  # Compute the refined matrix

  if((order - deriv) >= 1){
    design_all <- global_poly_helper(refined_x, p = order)
    design_all <- as.matrix(design_all[, 1:(order - deriv), drop = FALSE])
    for (i in 1:ncol(design_all)) {
      design_all[, i] <- (factorial(i + deriv - 1) / factorial(i - 1)) * design_all[, i]
    }
  }

  if (psd_iwp != 0) {
    B_refined <- local_poly_helper(knots = knots, refined_x = refined_x, p = (order-deriv))
    design_all <- cbind(B_refined, design_all)
  }else{
    B_refined <- matrix(0, nrow = length(refined_x), ncol = 0)
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

  # Add sigmaIWP for random effects
  if (psd_iwp != 0) {
    tmbdat$sigmaIWP <- psd_iwp / sqrt((pred_step ^ ((2 * order) - 1)) / (((2 * order) - 1) * (factorial(order - 1) ^ 2)))
  }

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

  # Generate posterior samples for coefficients
  samps_coef <- LaplacesDemon::rmvnp(n = M, mu = opt$par, Omega = as.matrix(prec_matrix))

  if(deriv > 0){
    if(ncol(B_refined) > 0){
      # take out the corresponding B elements and global poly elements
      samps_coef1 <- samps_coef[, 1:ncol(B_refined), drop = FALSE]
      samps_coef2 <- samps_coef[, (ncol(B_refined) + 1):ncol(samps_coef), drop = FALSE]
      # remove the first deriv columns
      samps_coef2 <- samps_coef2[, -c(1:deriv), drop = FALSE]
      # redefine samples_coef
      samps_coef <- cbind(samps_coef1, samps_coef2)
    }else{
      samps_coef <- samps_coef[, -c(1:deriv), drop = FALSE]
    }
  }

  # Compute posterior sample paths
  samps_fitted <- design_all %*% t(samps_coef)

  return(samps_fitted)
}





#' Perform Bayesian Model Averaging (BMA) Sampling
#'
#' This function performs Bayesian Model Averaging (BMA) for the i-th dataset by:
#' 1. Sampling PSD values based on posterior weights.
#' 2. Using the sampled PSD values to generate posterior samples.
#'
#' @param data_i A single dataset extracted from the \code{data_list} component of \code{fash_set_data}. Must be a list containing \code{y}, \code{x}, and \code{offset}.
#' @param posterior_weights A numeric vector representing the posterior weights of the PSD values for the i-th dataset.
#' @param psd_values A numeric vector of PSD values corresponding to the posterior weights.
#' @param refined_x A numeric vector of grid points where posterior sample paths are evaluated.
#' @param M An integer specifying the total number of posterior samples to generate.
#' @param Si A numeric vector representing the standard errors for the dataset (if applicable).
#' @param Omegai A numeric precision matrix for the dataset (if applicable).
#' @param num_basis An integer specifying the number of O-Spline basis functions.
#' @param betaprec A numeric value representing the precision of the fixed effects coefficients.
#' @param order An integer specifying the order of the Integrated Wiener Process (IWP) prior.
#' @param pred_step A numeric value specifying the prediction step size.
#' @param likelihood A character string specifying the likelihood function to use. Options are `gaussian` and `poisson`.
#' @param deriv An integer specifying the order of the derivative to compute.
#'
#' @return A list containing:
#'   \describe{
#'     \item{\code{posterior_samples}}{A numeric matrix of posterior sample paths, where rows correspond to grid points in \code{refined_x} and columns correspond to the M generated samples.}
#'     \item{\code{sampled_counts}}{A table of sampled PSD values and their corresponding frequencies.}
#'   }
#'
#' @examples
#' # Example usage
#' set.seed(1)
#' Y <- matrix(rpois(20, lambda = 5), nrow = 4, ncol = 5)
#' smooth_var <- matrix(runif(20), nrow = 4, ncol = 5)
#' offset <- 1
#' fash_data <- fash_set_data(Y, smooth_var, offset)
#' posterior_weights <- c(0.1, 0.4, 0.3, 0.2)
#' psd_values <- c(0.1, 0.5, 1.0, 2.0)
#' refined_x <- seq(0, 10, length.out = 50)
#' results <- fashr:::fash_bma_sampling(
#'   data_i = fash_data$data_list[[1]],
#'   posterior_weights = posterior_weights,
#'   psd_values = psd_values,
#'   refined_x = refined_x,
#'   M = 100,
#'   likelihood = "poisson"
#' )
#' print(results$posterior_samples)
#' print(results$sampled_counts)
#'
#' @keywords internal
#'
fash_bma_sampling <- function(data_i, posterior_weights, psd_values, refined_x, M, Si = NULL, Omegai = NULL,
                              num_basis = 30, betaprec = 1e-6, order = 1, pred_step = 1, likelihood, deriv = 0) {
  # Check that posterior_weights and psd_values match in length
  if (length(posterior_weights) != length(psd_values)) {
    stop("posterior_weights and psd_values must have the same length.")
  }

  # Normalize posterior_weights to ensure they sum to 1
  posterior_weights <- posterior_weights / sum(posterior_weights)

  # Sample PSD values (by index, so a single PSD value is never mistaken
  # for the sample() shorthand of 1:n) and count occurrences
  sampled_idx <- sample.int(length(psd_values), size = M, replace = TRUE, prob = posterior_weights)
  sampled_counts <- table(psd_values[sampled_idx])

  # Initialize matrix to store posterior samples
  n_points <- length(refined_x)
  posterior_samples <- matrix(0, nrow = n_points, ncol = M)

  # Loop through unique PSD values using sampled counts
  start_col <- 1
  for (psd in as.numeric(names(sampled_counts))) {
    # Number of samples for the current PSD
    num_samples <- sampled_counts[as.character(psd)]

    # Generate samples using fash_fit_once
    samples_for_psd <- fash_fit_once(
      data_i = data_i,
      refined_x = refined_x,
      M = num_samples,
      psd_iwp = psd,
      Si = Si,
      Omegai = Omegai,
      num_basis = num_basis,
      betaprec = betaprec,
      order = order,
      pred_step = pred_step,
      likelihood = likelihood,
      deriv = deriv
    )

    # Assign samples to the corresponding columns in the output matrix
    posterior_samples[, start_col:(start_col + num_samples - 1)] <- samples_for_psd
    start_col <- start_col + num_samples
  }

  return(list(posterior_samples = posterior_samples, sampled_counts = sampled_counts))

}


