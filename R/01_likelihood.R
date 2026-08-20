#' Compute the L Matrix for a Set of Datasets
#'
#' Computes the L matrix, where each row corresponds to a dataset and each column corresponds to a grid value of PSD (Predictive Standard Deviation).
#' Handles both Gaussian and Poisson likelihoods using helper functions.
#'
#' @param fash_data The output from \code{fash_set_data}, containing preprocessed datasets.
#' @param likelihood A character string specifying the likelihood function to use. Options are `gaussian` and `poisson`.
#' @param num_cores An integer specifying the number of cores to use for parallel processing.
#' @param grid A numeric vector representing the grid of PSD (Predictive Standard Deviation) values.
#' @param pred_step A numeric value specifying the prediction step size.
#' @param num_basis An integer specifying the number of O-Spline basis functions to use for the approximation.
#' @param betaprec A numeric value representing the precision of the fixed effects coefficients.
#' @param order An integer specifying the order of the Integrated Wiener Process (IWP) prior.
#' @param verbose A logical value. If \code{TRUE}, shows a progress bar when \code{num_cores = 1}.
#'
#' @return A numeric matrix where each row corresponds to a dataset and each column corresponds to a grid value of PSD.
#'
#' @examples
#' # Example usage
#' set.seed(1)
#' Y <- matrix(rpois(20, lambda = 5), nrow = 4, ncol = 5)
#' S <- matrix(1, nrow = 4, ncol = 5)
#' smooth_var <- matrix(runif(20), nrow = 4, ncol = 5)
#' offset <- 1
#' fash_data_poisson <- fash_set_data(Y = Y, smooth_var = smooth_var, offset = offset)
#' fash_data_gaussian <- fash_set_data(Y = Y, smooth_var = smooth_var, S = S)
#' grid <- seq(0.1, 2, length.out = 10)
#' L_matrix_gaussian <-
#'   fashr:::fash_L_compute(fash_data_gaussian, likelihood = "gaussian",
#'                          num_cores = 1, grid = grid, num_basis = 30,
#'                          betaprec = 1e-6, order = 2, verbose = TRUE)
#' L_matrix_poisson <-
#'   fashr:::fash_L_compute(fash_data_poisson, likelihood = "poisson",
#'                          num_cores = 1, grid = grid, num_basis = 30,
#'                          betaprec = 1e-6, order = 2, verbose = TRUE)
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom parallel mclapply
#'
#' @keywords internal
#'
fash_L_compute <- function(fash_data, likelihood = "gaussian", num_cores = 1, grid = seq(0, 2, length.out = 10), pred_step = 1, num_basis = 30, betaprec = 1e-6, order = 1, verbose = FALSE) {
  num_cores <- sanitize_num_cores(num_cores)
  num_datasets <- length(fash_data$data_list)
  datasets <- fash_data$data_list
  S_list <- fash_data$S
  Omega_list <- fash_data$Omega

  if (likelihood == "gaussian") {

    # Check if either S_list or Omega_list is not NULL
    if (is.null(S_list) && is.null(Omega_list)) {
      stop("Both S_list and Omega_list are NULL. At least one must be specified when likelihood = `gaussian`.")
    }


    if (num_cores == 1 && verbose) {
      # Sequential execution with progress bar
      pb <- utils::txtProgressBar(min = 0, max = num_datasets, style = 3)
      L_vecs <- lapply(1:num_datasets, function(i) {
        utils::setTxtProgressBar(pb, i)  # Update progress bar
        compute_L_gaussian_helper_seq(
          data_i = datasets[[i]],
          Si = S_list[[i]],
          Omegai = Omega_list[[i]],
          grid = grid,
          pred_step = pred_step,
          num_basis = num_basis,
          betaprec = betaprec,
          order = order
        )
      })
      close(pb)
    } else {
      # Parallel execution (no progress bar)
      L_vecs <- parallel::mclapply(1:num_datasets, function(i) {
        compute_L_gaussian_helper_seq(
          data_i = datasets[[i]],
          Si = S_list[[i]],
          Omegai = Omega_list[[i]],
          grid = grid,
          pred_step = pred_step,
          num_basis = num_basis,
          betaprec = betaprec,
          order = order
        )
      }, mc.cores = num_cores)
    }

    # Combine the L vectors into a matrix
    L_matrix <- do.call(rbind, L_vecs)
  }
  else if (likelihood == "poisson") {
    if (num_cores == 1 && verbose) {
      # Sequential execution with progress bar
      pb <- utils::txtProgressBar(min = 0, max = num_datasets, style = 3)
      L_vecs <- lapply(1:num_datasets, function(i) {
        utils::setTxtProgressBar(pb, i)  # Update progress bar
        compute_L_poisson_helper_seq(
          data_i = datasets[[i]],
          grid = grid,
          pred_step = pred_step,
          num_basis = num_basis,
          betaprec = betaprec,
          order = order
        )
      })
      close(pb)
    } else {
      # Parallel execution (no progress bar)
      L_vecs <- parallel::mclapply(1:num_datasets, function(i) {
        compute_L_poisson_helper_seq(
          data_i = datasets[[i]],
          grid = grid,
          pred_step = pred_step,
          num_basis = num_basis,
          betaprec = betaprec,
          order = order
        )
      }, mc.cores = num_cores)
    }

    # Combine the L vectors into a matrix
    L_matrix <- do.call(rbind, L_vecs)
  }
  else {
    stop("Unknown likelihood function. Currently, only 'gaussian' or 'poisson' is supported.")
  }

  return(L_matrix)
}





#' Compute Log-Likelihood for a Gaussian TMB Model
#'
#' Computes the log-likelihood for a Gaussian model using a TMB-based approach.
#' Handles both random-effects models and fixed-effects-only models, with DLL selection
#' based on whether standard errors (\code{S}) are specified.
#'
#' @param data_i A single dataset extracted from the \code{data_list} component of \code{fash_set_data}.
#'               Must be a list containing \code{y}, \code{x}, and \code{offset}.
#' @param Si A numeric vector representing the standard errors for the dataset.
#' @param Omegai A numeric precision matrix for the dataset.
#' @param psd_iwp A numeric value for the precision parameter of the Integrated Wiener Process. If 0, only fixed effects are used.
#' @param num_basis An integer specifying the number of O-Spline basis functions. Default is 30.
#' @param betaprec A numeric value representing the precision of the fixed effects coefficients.
#' @param order An integer specifying the order of the Integrated Wiener Process (IWP) prior.
#' @param pred_step A numeric value specifying the prediction step size.
#'
#' @return A numeric value representing the negative log-likelihood for the dataset.
#'
#' @examples
#' # Example usage
#' #' set.seed(1)
#' Y <- matrix(rnorm(20), nrow = 4, ncol = 5)
#' smooth_var <- matrix(runif(20), nrow = 4, ncol = 5)
#' offset <- 1
#' S <- c(0.5, 0.8, 1.2, 1.0, 0.9)
#' Omega <- diag(5)
#' data <- fash_set_data(Y, smooth_var, offset, S, Omega)
#' log_likelihood <- fashr:::compute_L_gaussian_helper(
#'   data$data_list[[1]], Si = S[1], Omegai = Omega, psd_iwp = 0.1
#' )
#'
#' @importFrom TMB MakeADFun
#' @keywords internal
#'
compute_L_gaussian_helper <- function(data_i, Si, Omegai, psd_iwp, num_basis = 30, betaprec = 1e-6, order = 1, pred_step = 1) {
  # Create the tmbdat object using existing helper function
  tmbdat <- fash_set_tmbdat(data_i, Si, Omegai, num_basis = num_basis, betaprec = betaprec, order = order)

  if (psd_iwp != 0) {
    # Add sigmaIWP for random effects
    tmbdat$sigmaIWP <- psd_iwp / sqrt((pred_step ^ ((2 * order) - 1)) / (((2 * order) - 1) * (factorial(order - 1) ^ 2)))

    # Define initial parameter values
    tmbparams <- list(
      W = rep(0, ncol(tmbdat$X) + ncol(tmbdat$B))
    )

    # Choose the appropriate DLL based on whether S is specified
    DLL <- if (!is.null(Si)) {
      "Gaussian_ind"
    } else {
      "Gaussian_dep"
    }

    # Create the TMB model
    ff <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      random = "W",
      DLL = DLL,
      silent = TRUE
    )
  } else {
    # Fixed effects only (no random effects)
    tmbparams <- list(
      W = rep(0, ncol(tmbdat$X))
    )

    # Choose the appropriate DLL based on whether S is specified
    DLL <- if (!is.null(Si)) {
      "Gaussian_ind_fixed"
    } else {
      "Gaussian_dep_fixed"
    }

    # Create the TMB model
    ff <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      random = "W",
      DLL = DLL,
      silent = TRUE
    )
  }

  # Compute and return the negative log-likelihood
  return(-ff$fn())
}





#' Compute Log-Likelihoods for All Grid Values of PSD (Predictive Standard Deviation)
#'
#' Computes the log-likelihood for a single dataset across all values of \code{psd_iwp} specified in the grid.
#' This function iteratively calls \code{compute_L_gaussian_helper} for each grid value.
#'
#' @param data_i A single dataset extracted from the \code{data_list} component of \code{fash_set_data}.
#'               Must be a list containing \code{y}, \code{x}, and \code{offset}.
#' @param Si A numeric vector representing the standard errors for the dataset.
#' @param Omegai A numeric precision matrix for the dataset.
#' @param grid A numeric vector representing the grid of PSD (Predictive Standard Deviation) values.
#' @param num_basis An integer specifying the number of O-Spline basis functions.
#' @param betaprec A numeric value representing the precision of the fixed effects coefficients.
#' @param order An integer specifying the order of the Integrated Wiener Process (IWP) prior.
#' @param pred_step A numeric value specifying the prediction step size.
#'
#' @return A numeric vector of log-likelihood values, one for each grid value.
#'
#' @examples
#' # Example usage
#' set.seed(1)
#' Y <- matrix(rnorm(20), nrow = 4, ncol = 5)
#' smooth_var <- matrix(runif(20), nrow = 4, ncol = 5)
#' offset <- 1
#' S <- c(0.5, 0.8, 1.2, 1.0, 0.9)
#' Omega <- diag(5)
#' data <- fash_set_data(Y, smooth_var, offset, S, Omega)
#' grid <- seq(0.1, 2, length.out = 10)
#' likelihoods <- fashr:::compute_L_gaussian_helper_seq(
#'   data$data_list[[1]], Si = S[[1]], Omegai = Omega, grid = grid
#' )
#'
#' @keywords internal
#'
compute_L_gaussian_helper_seq <- function(data_i, Si, Omegai, grid, num_basis = 30, betaprec = 1e-6, order = 1, pred_step = 1) {
  # Initialize vector to store log-likelihoods
  log_likelihoods <- numeric(length(grid))

  # Loop over grid values
  for (g in seq_along(grid)) {
    log_likelihoods[g] <- compute_L_gaussian_helper(
      data_i = data_i,
      Si = Si,
      Omegai = Omegai,
      psd_iwp = grid[g],
      num_basis = num_basis,
      betaprec = betaprec,
      order = order,
      pred_step = pred_step
    )
  }

  return(log_likelihoods)
}





#' Compute Log-Likelihoods for All Grid Values of PSD (Poisson)
#'
#' Computes the log-likelihood for a single dataset across all values of \code{psd_iwp} specified in the grid.
#' This function iteratively calls \code{compute_L_poisson_helper} for each grid value.
#'
#' @param data_i A single dataset extracted from the \code{data_list} component of \code{fash_set_data}.
#'               Must be a list containing \code{y}, \code{x}, and \code{offset}.
#' @param grid A numeric vector representing the grid of PSD (Predictive Standard Deviation) values.
#' @param num_basis An integer specifying the number of O-Spline basis functions.
#' @param betaprec A numeric value representing the precision of the fixed effects coefficients (`beta`).
#' @param order An integer specifying the order of the Integrated Wiener Process (IWP) prior.
#' @param pred_step A numeric value specifying the prediction step size.
#'
#' @return A numeric vector of log-likelihood values, one for each grid value.
#'
#' @examples
#' # Example usage
#' set.seed(1)
#' Y <- matrix(rpois(20, lambda = 5), nrow = 4, ncol = 5)
#' smooth_var <- matrix(runif(20), nrow = 4, ncol = 5)
#' offset <- 1
#' fash_data <- fash_set_data(Y, smooth_var, offset)
#' grid <- seq(0.1, 2, length.out = 10)
#' likelihoods <- fashr:::compute_L_poisson_helper_seq(
#'   data_i = fash_data$data_list[[1]], grid = grid
#' )
#'
#' @keywords internal
#'
compute_L_poisson_helper_seq <- function(data_i, grid, num_basis = 30, betaprec = 1e-6, order = 1, pred_step = 1) {
  # Initialize vector to store log-likelihoods
  log_likelihoods <- numeric(length(grid))

  # Loop over grid values
  for (g in seq_along(grid)) {
    log_likelihoods[g] <- compute_L_poisson_helper(
      data_i = data_i,
      psd_iwp = grid[g],
      num_basis = num_basis,
      betaprec = betaprec,
      order = order,
      pred_step = pred_step
    )
  }

  return(log_likelihoods)
}








#' Compute Log-Likelihood for a Poisson TMB Model
#'
#' Computes the log-likelihood for a Poisson model using a TMB-based approach.
#' Handles both random-effects models and fixed-effects-only models.
#'
#' @param data_i A single dataset extracted from the \code{data_list} component of \code{fash_set_data}.
#'               Must be a list containing \code{y}, \code{x}, and \code{offset}.
#' @param psd_iwp A numeric value for the precision parameter of the Integrated Wiener Process. If 0, only fixed effects are used.
#' @param num_basis An integer specifying the number of O-Spline basis functions.
#' @param betaprec A numeric value representing the precision of the fixed effects coefficients (`beta`).
#' @param order An integer specifying the order of the Integrated Wiener Process (IWP) prior.
#' @param pred_step A numeric value specifying the prediction step size.
#'
#' @return A numeric value representing the negative log-likelihood for the dataset.
#'
#' @examples
#' # Example usage
#' set.seed(1)
#' Y <- matrix(rpois(20, lambda = 5), nrow = 4, ncol = 5)
#' smooth_var <- matrix(runif(20), nrow = 4, ncol = 5)
#' offset <- 1
#' fash_data <- fash_set_data(Y, smooth_var, offset)
#' log_likelihood <- fashr:::compute_L_poisson_helper(
#'   data_i = fash_data$data_list[[1]], psd_iwp = 0.1
#' )
#'
#' @importFrom TMB MakeADFun
#'
#' @keywords internal
#'
compute_L_poisson_helper <- function(data_i, psd_iwp, num_basis = 30, betaprec = 1e-6, order = 1, pred_step = 1) {
  # Create the tmbdat object using existing helper function
  tmbdat <- fash_set_tmbdat(data_i, num_basis = num_basis, betaprec = betaprec, order = order)

  if (psd_iwp != 0) {
    # Add sigmaIWP for random effects
    tmbdat$sigmaIWP <- psd_iwp / sqrt((pred_step ^ ((2 * order) - 1)) / (((2 * order) - 1) * (factorial(order - 1) ^ 2)))

    # Define initial parameter values
    tmbparams <- list(
      W = rep(0, ncol(tmbdat$X) + ncol(tmbdat$B))
    )

    DLL <- "Poisson_ind"

    # Create the TMB model
    ff <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      random = "W",
      DLL = DLL,
      silent = TRUE
    )
  } else {
    # Fixed effects only (no random effects)
    tmbparams <- list(
      W = rep(0, ncol(tmbdat$X))
    )

    DLL <- "Poisson_ind_fixed"

    # Create the TMB model
    ff <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      random = "W",
      DLL = DLL,
      silent = TRUE
    )
  }

  # Compute and return the negative log-likelihood
  return(-ff$fn())
}









