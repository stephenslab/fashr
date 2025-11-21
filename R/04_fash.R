#' Perform Full FASH Analysis
#'
#' This function performs the full FASH pipeline, including data setup, likelihood computation,
#' empirical Bayes estimation, and outputs a structured \code{fash} object.
#'
#' @param Y Either a numeric matrix of response variables or a character string specifying the column name in \code{data_list} for response variables.
#' @param smooth_var A numeric matrix, vector, or a character string specifying the column name in \code{data_list} for smoothing variables.
#' @param offset A numeric matrix, vector, scalar, or a character string specifying the column name in \code{data_list} for offset variables.
#' @param S A numeric matrix, vector, scalar, or list representing the standard errors of \code{Y}. Or a character string specifying the column name in \code{data_list} for SD.
#' @param Omega Either a list of precision matrices (one for each dataset) or a single precision matrix (shared across all datasets).
#' @param data_list A list of data frames, where each data frame corresponds to a single dataset.
#' @param grid A numeric vector representing the grid of PSD (Predictive Standard Deviation) values.
#' @param likelihood A character string specifying the likelihood function to use. Options are `gaussian` and `poisson`.
#' @param num_basis An integer specifying the number of O-Spline basis functions.
#' @param betaprec A numeric value representing the precision of the fixed effects coefficients.
#' @param order An integer specifying the order of the Integrated Wiener Process (IWP) prior.
#' @param pred_step A numeric value specifying the prediction step size.
#' @param penalty A numeric value representing the lambda value for the Dirichlet prior.
#' @param num_cores An integer specifying the number of cores to use for parallel processing.
#' @param verbose A logical value. If \code{TRUE}, shows progress messages and timing for each step.
#'
#' @return A \code{fash} object containing:
#'   \describe{
#'     \item{\code{prior_weights}}{Estimated prior weights for PSD values.}
#'     \item{\code{posterior_weights}}{Posterior weight matrix of PSD values.}
#'     \item{\code{psd_grid}}{PSD grid values.}
#'     \item{\code{lfdr}}{Local false discovery rate for each dataset.}
#'     \item{\code{settings}}{A list of settings used in the FASH pipeline.}
#'     \item{\code{fash_data}}{A structured list of data components.}
#'     \item{\code{L_matrix}}{Likelihood matrix used in the FASH pipeline.}
#'     \item{\code{eb_result}}{Empirical Bayes estimation results.}
#'   }
#'
#' @examples
#' set.seed(1)
#' data_list <- list(
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0),
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0)
#' )
#' grid <- seq(0, 2, length.out = 10)
#' result <- fash(data_list = data_list, Y = "y", smooth_var = "x", offset = "offset", grid = grid, likelihood = "poisson", verbose = TRUE)
#' print(result)
#'
#' @importFrom graphics abline
#' @importFrom graphics legend
#'
#' @export
fash <- function(Y, smooth_var, offset = 0, S = NULL,
                 Omega = NULL, data_list = NULL,
                 grid = seq(0, 1, length.out = 25),
                 likelihood = c("gaussian","poisson"),
                 num_basis = 30, betaprec = 1e-6, order = 1, pred_step = 1,
                 penalty = 1, num_cores = 1, verbose = TRUE) {

  # Check if order is a positive integer
  if (!is.numeric(order) || length(order) != 1 || order <= 0 || order != floor(order)) {
    stop("Order must be a positive integer.")
  }

  # If order is larger than 4, give a warning about slower computation
  if (order > 4) {
    # print a message
    if (verbose) {
      cat(sprintf("Note: Order = %d. Large choice of order may lead to slower computation times.\n", order))
    }
  }


  # Check if 0 is included in the grid, if not add it and produce a warning
  if (!0 %in% grid) {
    warning("0 is not included in the grid, adding it to the grid.")
    grid <- c(0, grid)
  }

  # If likelihood is "gaussian", ensure either S or Omega is provided
  likelihood <- match.arg(likelihood)
  if (likelihood == "gaussian") {
    if (is.null(S) && is.null(Omega)) {
      warning("For Gaussian likelihood, either S or Omega must be provided. Defaulting to S = 1 for all datasets.")
      S <- 1
    }
    if (!is.null(S) && !is.null(Omega)) {
      warning("Both S and Omega are provided. Using S for standard errors.")
    }
  }

  # Helper function for timing and verbose output
  timing_message <- function(step_name, code_block) {
    start_time <- Sys.time()
    if (verbose) cat(sprintf("Starting %s...\n", step_name))
    result <- code_block()
    elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")

    # Format elapsed time dynamically
    if (elapsed_time < 60) {
      time_message <- sprintf("%.2f seconds", elapsed_time)
    } else if (elapsed_time < 3600) {
      time_message <- sprintf("%.2f minutes", elapsed_time / 60)
    } else {
      time_message <- sprintf("%.2f hours", elapsed_time / 3600)
    }

    if (verbose) cat(sprintf("Completed %s in %s.\n", step_name, time_message))
    return(result)
  }

  # Step 1: Data setup
  fash_data <- timing_message("data setup", function() {
    fash_set_data(data_list = data_list, Y = Y, smooth_var = smooth_var, offset = offset, S = S, Omega = Omega)
  })

  # Step 2: Likelihood computation
  L_matrix <- timing_message("likelihood computation", function() {
    fash_L_compute(fash_data, likelihood = likelihood, num_cores = num_cores, grid = grid,
                   num_basis = num_basis, betaprec = betaprec, order = order, pred_step = pred_step,
                   verbose = verbose)
  })

  # Step 3: Empirical Bayes estimation
  eb_result <- timing_message("empirical Bayes estimation", function() {
    fash_eb_est(L_matrix, grid = grid, penalty = penalty)
  })
  # Defining rownames of posterior weights
  rownames(eb_result$posterior_weight) <- names(fash_data$data_list)  # Ensure dataset names are set

  # Add dataset names if missing
  if (is.null(rownames(eb_result$posterior_weight))) {
    rownames(eb_result$posterior_weight) <- paste0("Dataset_", seq_len(nrow(eb_result$posterior_weight)))
  }


  # Step 4: Compute additional metrics
  # if psd_value zero is included:
  if(0 %in% eb_result$prior_weight$psd){
    lfdr <- eb_result$posterior_weight[, which(eb_result$prior_weight$psd == 0)]
  }else{
    lfdr <- rep(0, nrow(eb_result$posterior_weight))
  }

  # Step 5: Create and return the fash object
  result <- structure(
    list(
      prior_weights = eb_result$prior_weight,
      posterior_weights = eb_result$posterior_weight,
      psd_grid = grid,
      lfdr = lfdr,
      settings = list(
        num_basis = num_basis,
        betaprec = betaprec,
        order = order,
        pred_step = pred_step,
        likelihood = likelihood,
        penalty = penalty
      ),
      fash_data = fash_data,
      L_matrix = L_matrix,
      eb_result = eb_result
    ),
    class = "fash"
  )

  if (verbose) cat("fash object created successfully.\n")

  return(result)
}







#' Perform False Discovery Rate (FDR) Control
#'
#' This function performs hypothesis testing by controlling the False Discovery Rate (FDR) based on the
#' local false discovery rate (lfdr) stored in the \code{fash} object.
#'
#' @param fash_obj A \code{fash} object containing the results of the FASH pipeline,
#'   including the vector \code{lfdr}.
#' @param alpha A numeric value specifying the FDR threshold (significance level)
#'   used to declare discoveries. Default is \code{0.05}.
#' @param plot A logical value. If \code{TRUE}, plots the cumulative FDR values
#'   against the rank of units sorted by lfdr, with a horizontal line at the
#'   \code{alpha} level. Default is \code{FALSE}.
#' @param sort A logical value. If \code{TRUE}, returns the FDR results sorted by increasing lfdr.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{fdr_results}{A data frame with one row per unit,
#'     containing:
#'     \itemize{
#'       \item \code{rank}: Rank of the unit after sorting by lfdr (1 = smallest lfdr).
#'       \item \code{index}: Original index of the unit in \code{lfdr}.
#'       \item \code{lfdr}: The sorted local false discovery rate values.
#'       \item \code{FDR}: The cumulative FDR at each rank, computed as
#'         the running mean of the sorted lfdr values.
#'     }}
#'   \item{message}{A character string summarizing how many units are significant
#'     at the chosen \code{alpha} level and the total number of units tested.}
#'   \item{significant_units}{A vector of indices (or names, if \code{lfdr} is named)
#'     corresponding to units with cumulative FDR less than or equal to \code{alpha}.}
#' }
#'
#' @examples
#' # Example usage
#' set.seed(1)
#' data_list <- list(
#'   data.frame(y = rnorm(n = 5, sd = 0.5), x = 1:5, offset = 0),
#'   data.frame(y = rnorm(n = 5, sd = 0.8), x = 1:5, offset = 0),
#'   data.frame(y = rnorm(n = 5, sd = 0.6), x = 1:5, offset = 0),
#'   data.frame(y = rnorm(n = 5, sd = 0.7), x = 1:5, offset = 0)
#' )
#' S <- list(rep(0.5, 5), rep(0.8, 5), rep(0.6, 5), rep(0.7, 5))
#' grid <- seq(0, 2, length.out = 10)
#' fash_obj <- fash(data_list = data_list, Y = "y", smooth_var = "x", offset = "offset", S = S, grid = grid, likelihood = "gaussian", verbose = TRUE)
#' fdr_control(fash_obj, alpha = 0.05, plot = TRUE)
#'
#' @export
fdr_control <- function(fash_obj, alpha = 0.05, plot = FALSE, sort = FALSE) {
  # Extract lfdr
  lfdr <- fash_obj$lfdr
  if (is.null(lfdr)) {
    stop("The `fash` object does not contain local false discovery rates (lfdr).")
  }

  # Sort lfdr
  n <- length(lfdr)
  lfdr_sorted <- sort(lfdr, index.return = TRUE)
  cumulative_lfdr <- cumsum(lfdr_sorted$x) / seq_len(n)

  # Identify significant datasets
  significant <- which(cumulative_lfdr <= alpha)
  significant_count <- length(significant)

  # Prepare results (sorted table)
  fdr_results <- data.frame(
    rank = seq_len(n),
    index = lfdr_sorted$ix,
    lfdr = lfdr_sorted$x,
    FDR  = cumulative_lfdr
  )

  # If sort is FALSE, reorder fdr_results to original order
  if (!sort) {
    fdr_results <- fdr_results[order(fdr_results$index), ]
  }

  # Prepare significant units (in original indices or names)
  significant_units <- lfdr_sorted$ix[significant]
  if (!is.null(names(lfdr))) {
    significant_units <- names(lfdr)[significant_units]
  }

  # Display message
  message <- sprintf(
    "%d datasets are significant at alpha level %.2f. Total datasets tested: %d.",
    significant_count, alpha, n
  )
  cat(message, "\n")

  # Plot
  if (plot) {
    plot(
      1:n, cumulative_lfdr, type = "b", pch = 19, col = "blue",
      xlab = "Dataset Rank (Sorted by LFDR)", ylab = "Cumulative FDR",
      main = sprintf("FDR Control with Alpha = %.2f", alpha)
    )
    abline(h = alpha, col = "red", lty = 2)
    legend(
      "topright",
      legend = c("FDR Values", "Alpha Level"),
      col = c("blue", "red"), lty = c(1, 2), pch = c(19, NA)
    )
  }

  # Return results
  list(
    fdr_results      = fdr_results,
    message          = message,
    significant_units = significant_units
  )
}


#' Plot Method for fash Objects
#'
#' Generates a plot for a \code{fash} object, providing either a heatmap of posterior weights
#' or a structure plot summarizing component contributions across datasets.
#'
#' @param x A \code{fash} object containing the results of the FASH pipeline.
#' @param plot_type A character string specifying the type of plot to generate.
#'   One of:
#'   - \code{"heatmap"}: Bubble/heatmap plot of posterior weights (default).
#'   - \code{"structure"}: Structure plot of mixture components.
#'   - \code{"function"}: Plot fitted effect function for a selected unit.
#' @param ordering A character string specifying the method for reordering datasets in the structure plot.
#'   Only used if \code{plot_type = "structure"}.
#'
#'   - \code{"mean"}: Reorder by the mean of the posterior PSD.
#'   - \code{"lfdr"}: Reorder by the local false discovery rate.
#'   - \code{NULL}: No reordering (default).
#'
#' @param discrete A logical value. If \code{TRUE}, treats PSD values as discrete categories with distinct colors
#'                 in the structure plot. Ignored if \code{plot_type = "heatmap"} or \code{"function"}.
#' @param ... Additional arguments passed to \code{plot_heatmap}, \code{fash_structure_plot} or \code{plot_function},
#'
#' @return A plot object (typically a \code{ggplot}).
#'
#' @examples
#' set.seed(1)
#' data_list <- list(
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0),
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0)
#' )
#' grid <- seq(0, 2, length.out = 10)
#' fash_obj <- fash(data_list = data_list, Y = "y", smooth_var = "x", offset = "offset", grid = grid, likelihood = "poisson", verbose = TRUE)
#'
#' # Heatmap plot
#' plot(fash_obj)
#'
#' # Structure plot
#' plot(fash_obj, plot_type = "structure", ordering = "mean", discrete = TRUE)
#'
#' @export
plot.fash <- function(x,
                      plot_type = c("heatmap", "structure", "function"),
                      ordering = NULL,
                      discrete = FALSE,
                      selected_unit = NULL,
                      ...) {
  # Match the plot type
  plot_type <- match.arg(plot_type)

  # Validate input
  if (!inherits(x, "fash")) {
    stop("Input must be a `fash` object.")
  }

  if (plot_type == "heatmap") {
    return(
      plot_heatmap(
        object = x,
        ...
      )
    )
  }

  else if (plot_type == "structure") {
    return(
      fash_structure_plot(
        eb_output = list(
          posterior_weight = x$posterior_weights,
          prior_weight = x$prior_weights
        ),
        ordering = ordering,
        discrete = discrete,
        ...
      )
    )
  }

  else if (plot_type == "function") {
    if (is.null(selected_unit)) {
      stop("Please provide 'selected_unit' for function plot.")
    }
    return(
      plot_function(
        fash_obj = x,
        selected_unit = selected_unit,
        ...
      )
    )
  }

  else {
    stop("Invalid plot_type. Must be either 'heatmap', 'structure', or 'function'.")
  }
}


#' Plot Fitted Effect Function for a Selected Unit
#'
#' This internal helper function generates a plot for a single unit
#' in a fitted \code{fash} object. It overlays the observed data with the
#' smoothed posterior mean function and its uncertainty band.
#'
#' @param fash_obj A \code{fash} object produced by the \code{fash()} pipeline.
#' @param selected_unit An integer specifying which unit (dataset index)
#'   to visualize.
#' @param smooth_var Optional numeric vector giving the values of the
#'   smoothing variable at which the effect function should be evaluated.
#'   If \code{NULL}, 100 equally spaced points spanning the observed
#'   \code{x}-range of the selected unit are generated automatically.
#' @param ... Additional arguments passed to the underlying \code{plot()} call.
#'
#' @details
#' The function extracts the observed data for the selected unit and evaluates
#' the fitted effect function using \code{predict()}. The resulting mean
#' function and its 95\% credible band are plotted together with the raw data.
#'
#' @return
#' This function is called for generating a plot
#' and does not return a value.
#'
#' @keywords internal
plot_function <- function(fash_obj, selected_unit, smooth_var = NULL, ...) {
  # Extract dataset
  dataset <- fash_obj$fash_data$data_list[[selected_unit]]

  # If smooth_var is NULL, generate a sequence spanning the observed x-range
  if (is.null(smooth_var)) {
    smooth_var <- seq(min(dataset$x), max(dataset$x), length.out = 100)
  }

  # Get fitted effect
  fitted_beta_new <- predict(
    fash_obj,
    index = selected_unit,
    smooth_var = smooth_var
  )

  # Plot observed data
  plot(
    dataset$x, dataset$y,
    type = "p", col = "black", lwd = 2,
    xlab = "Condition", ylab = "Effect Size",
    main = paste("Unit", selected_unit),
    ...
  )

  # Add posterior mean curve
  lines(
    fitted_beta_new$x, fitted_beta_new$mean,
    col = "red", lwd = 2
  )

  # Add credible interval band
  polygon(
    c(fitted_beta_new$x, rev(fitted_beta_new$x)),
    c(fitted_beta_new$lower, rev(fitted_beta_new$upper)),
    col = rgb(1, 0, 0, 0.2), border = NA
  )
}



#' Predict Method for fash Objects
#'
#' Generates posterior predictions for a specific dataset from a \code{fash} object using Bayesian Model Averaging.
#'
#' @param object A \code{fash} object containing the results of the FASH pipeline.
#' @param index An integer specifying the dataset index to predict.
#' @param smooth_var A numeric vector specifying refined x values for prediction. If \code{NULL}, uses the x values from the model fit.
#' @param only.samples A logical value. If \code{TRUE}, returns posterior samples. If \code{FALSE}, summarizes the samples into mean and 95 percent confidence intervals.
#' @param M An integer specifying the number of posterior samples to generate.
#' @param deriv An integer specifying the order of the derivative to compute.
#' @param ... Additional arguments (not used).
#'
#' @return If \code{only.samples = TRUE}, a matrix of posterior samples where rows correspond to \code{smooth_var} and columns correspond to posterior draws.
#' If \code{only.samples = FALSE}, a data frame summarizing posterior predictions with columns:
#'
#'   - \code{x}: The refined x values.
#'
#'   - \code{mean}: The posterior mean.
#'
#'   - \code{lower}: The lower bound of the 95 percent interval.
#'
#'   - \code{upper}: The upper bound of the 95 percent interval.
#'
#'   - \code{median}: The posterior median.
#'
#' @examples
#'
#' set.seed(1)
#'
#' # Example 1: Predict for a specific dataset with summarized results
#' data_list <- list(
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0),
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0)
#' )
#' S <- list(rep(0.5, 5), rep(0.8, 5))
#' Omega <- list(diag(5), diag(5))
#' grid <- seq(0, 2, length.out = 10)
#' fash_obj <- fash(data_list = data_list, Y = "y", smooth_var = "x", offset = "offset", S = S, Omega = Omega, grid = grid, likelihood = "poisson", verbose = TRUE)
#' predict(fash_obj, index = 1, smooth_var = seq(1, 5, length.out = 50), only.samples = FALSE)
#'
#' # Example 2: Generate posterior samples
#' samples <- predict(fash_obj, index = 1, smooth_var = seq(1, 5, length.out = 50), only.samples = TRUE)
#' dim(samples)  # Rows: refined_x, Columns: posterior samples
#'
#' # Example 3: Use original x values for prediction
#' summary <- predict(fash_obj, index = 1, smooth_var = NULL, only.samples = FALSE)
#' head(summary)
#'
#' # Example 4: Increase number of posterior samples
#' samples <- predict(fash_obj, index = 1, smooth_var = seq(1, 5, length.out = 50), only.samples = TRUE, M = 5000)
#' summary <- predict(fash_obj, index = 1, smooth_var = seq(1, 5, length.out = 50), only.samples = FALSE, M = 5000)
#'
#' @importFrom stats predict
#' @importFrom stats median
#' @importFrom stats quantile
#'
#' @method predict fash
#'
#' @export
#'
predict.fash <- function (object, index = 1, smooth_var = NULL, only.samples = FALSE, M = 3000, deriv = 0, ...) {
  # Validate input
  if (!inherits(object, "fash")) {
    stop("Input must be a `fash` object.")
  }

  if(is.numeric(index)){
    if (index < 1 || index > length(object$posterior_weights)) {
      stop("Index is out of range for the datasets in the `fash` object.")
    }
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

  # Generate posterior samples using BMA
  posterior_samples <- fash_bma_sampling(
    data_i = data_i,
    posterior_weights = posterior_weights,
    psd_values = psd_values,
    refined_x = refined_x,
    M = M,
    Si = Si,
    Omegai = Omegai,
    num_basis = settings$num_basis,
    betaprec = settings$betaprec,
    order = settings$order,
    pred_step = settings$pred_step,
    likelihood = settings$likelihood,
    deriv = deriv
  )$posterior_samples

  # Return samples if only.samples = TRUE
  if (only.samples) {
    return(posterior_samples)
  }

  # Summarize samples: mean and 95% intervals
  posterior_summary <- data.frame(
    x = refined_x,
    mean = rowMeans(posterior_samples),
    median = apply(posterior_samples, 1, median),
    lower = apply(posterior_samples, 1, function(x) quantile(x, 0.025)),
    upper = apply(posterior_samples, 1, function(x) quantile(x, 0.975))
  )

  return(posterior_summary)
}





#' Print Method for fash Objects
#'
#' Displays a summary of the fitted \code{fash} object, including the number of datasets,
#' type of likelihood used, the number of PSD grid values, and the order of the Integrated Wiener Process (IWP).
#'
#' @param x A \code{fash} object.
#' @param ... Additional arguments (not used).
#'
#' @examples
#' set.seed(1)
#' data_list <- list(
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0),
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0)
#' )
#' grid <- seq(0, 2, length.out = 10)
#' fash_obj <- fash(data_list = data_list, Y = "y", smooth_var = "x", offset = "offset", grid = grid, likelihood = "poisson", verbose = TRUE)
#' print(fash_obj)
#'
#' @export
print.fash <- function(x, ...) {
  # Validate input
  if (!inherits(x, "fash")) {
    stop("Input must be a `fash` object.")
  }

  # Extract relevant information
  n_datasets <- length(x$fash_data$data_list)
  likelihood <- x$settings$likelihood
  n_grid_initial <- length(x$psd_grid)
  nontrivial_grid_values <- sum(x$prior_weights$prior_weight > 0)
  iwp_order <- x$settings$order

  # Display summary
  cat("Fitted fash Object\n")
  cat("-------------------\n")
  cat(sprintf("Number of datasets: %d\n", n_datasets))
  cat(sprintf("Likelihood: %s\n", likelihood))
  cat(sprintf("Number of PSD grid values: %d (initial), %d (non-trivial)\n", n_grid_initial, nontrivial_grid_values))
  cat(sprintf("Order of Integrated Wiener Process (IWP): %d\n", iwp_order))
}




#' Perform Functional Hypothesis Testing on Posterior Samples
#'
#' This function applies a user-specified functional to posterior samples from a \code{fash} object, calculates the
#' local false sign rate (LFSR) for each dataset, and returns a ranked data frame. The computation can be
#' parallelized if \code{num_cores > 1}.
#'
#' @param functional A function applied to each posterior sample to extract a scalar statistic.
#' @param lfsr_cal A function used to compute the local false sign rate (lfsr).
#' @param fash A \code{fash} object.
#' @param indices A numeric vector specifying the dataset indices to evaluate.
#' @param smooth_var A numeric vector specifying refined x values for prediction.
#' @param num_cores An integer specifying the number of cores to use for parallel processing.
#'
#' @return A data frame containing:
#'
#' \describe{
#'   \item{indices}{The dataset indices corresponding to \code{indices}.}
#'   \item{lfsr}{The computed local false sign rate (LFSR) for each dataset.}
#'   \item{cfsr}{The cumulative false sign rate (CFSR), calculated as the cumulative mean of \code{lfsr}.}
#' }
#'
#'
#' @examples
#' set.seed(1)
#'
#' # Define a functional (e.g., mean of posterior samples)
#' functional_example <- function(x) { mean(x) }
#'
#' # Example fash object (assuming it has been fitted)
#' data_list <- list(
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0),
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0)
#' )
#' grid <- seq(0, 2, length.out = 10)
#' fash_obj <- fash(data_list = data_list, Y = "y", smooth_var = "x", grid = grid, likelihood = "poisson", verbose = TRUE)
#'
#' # Perform functional hypothesis testing with parallel execution
#' result <- testing_functional(functional = functional_example, fash = fash_obj, indices = 1:2, num_cores = 2)
#' print(result)
#'
#'
#' @importFrom parallel mclapply
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @export
#'
testing_functional <- function(functional,
                               lfsr_cal = function(x) { min(mean(x <= 0), mean(x >= 0)) },
                               fash, indices,
                               smooth_var = NULL,
                               num_cores = 1) {
  # Define the function to be run for each index
  compute_lfsr <- function(index) {
    sample_index <- predict(fash, index = index, smooth_var = smooth_var, only.samples = TRUE)
    result <- apply(sample_index, 2, functional)
    lfsr <- lfsr_cal(result)
    return(c(index, lfsr))
  }

  # Parallel or sequential execution
  if (num_cores > 1) {
    results_list <- parallel::mclapply(indices, compute_lfsr, mc.cores = num_cores)
  } else {
    # Sequential execution with progress bar
    lfsr_vec <- NULL
    pb <- utils::txtProgressBar(min = 0, max = length(indices), style = 3)
    results_list <- list()
    for (i in seq_along(indices)) {
      utils::setTxtProgressBar(pb, i)
      results_list[[i]] <- compute_lfsr(indices[i])
    }
    close(pb)
  }

  # Convert results to a data frame
  results_mat <- do.call(rbind, results_list)

  result_df <- data.frame(
    indices = results_mat[, 1],
    lfsr    = results_mat[, 2]
  )
  result_df <- result_df[order(result_df$lfsr), ]

  result_df$cfsr <- cumsum(result_df$lfsr) / seq_len(nrow(result_df))

  return(result_df)
}





#' Structure Plot for Posterior Weights
#'
#' This function takes the output of \code{fash_eb_est} and generates a structure plot
#' visualizing the posterior weights for all datasets. It can display PSD values
#' as either continuous or discrete variables and optionally reorder datasets.
#'
#' @param eb_output A list output from \code{fash_eb_est}, containing:
#'   \describe{
#'     \item{posterior_weight}{A numeric matrix of posterior weights (datasets as rows, PSD as columns).}
#'     \item{prior_weight}{A data frame of prior weights (not used in this plot).}
#'   }
#' @param discrete A logical value. If \code{TRUE}, treats PSD values as discrete categories with distinct colors.
#'                 If \code{FALSE}, treats PSD values as a continuous variable with a gradient.
#' @param ordering A character string specifying the method for reordering datasets. Options are:
#'   \describe{
#'     \item{NULL}{No reordering (default).}
#'     \item{`mean`}{Reorder by the mean of the posterior PSD.}
#'     \item{`median`}{Reorder by the median of the posterior PSD.}
#'     \item{`lfdr`}{Reorder by the local false discovery rate (posterior probability of PSD = 0).}
#'   }
#' @param selected_indices A numeric vector specifying the indices of datasets to display. If \code{NULL}, displays all datasets.
#' @return A ggplot object representing the structure plot.
#'
#' @examples
#' # Example usage
#' set.seed(1)
#' grid <- seq(0.1, 2, length.out = 5)
#' L_matrix <- matrix(rnorm(20), nrow = 4, ncol = 5)
#' eb_output <- fash_eb_est(L_matrix, penalty = 2, grid = grid)
#' plot_cont <- fash_structure_plot(eb_output, discrete = FALSE, ordering = "mean")
#' plot_disc <- fash_structure_plot(eb_output, discrete = TRUE, ordering = "median")
#' print(plot_cont)
#' print(plot_disc)
#'
#' @importFrom ggplot2 ggplot aes geom_bar labs scale_fill_brewer
#' @importFrom ggplot2 scale_fill_gradient coord_flip theme_minimal
#' @importFrom ggplot2 theme element_blank element_rect
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#'
#' @export
#'
fash_structure_plot <- function (eb_output, discrete = FALSE,
                                 ordering = NULL,
                                 selected_indices = NULL) {

  # Select indices if specified
  if (!is.null(selected_indices)) {
    eb_output_selected <- eb_output
    eb_output_selected$posterior_weight <- eb_output$posterior_weight[selected_indices, , drop = FALSE]
  } else {
    eb_output_selected <- eb_output
  }

  # Extract posterior weights matrix
  posterior_weights_matrix <- eb_output_selected$posterior_weight

  # Reorder datasets if ordering is specified
  if (!is.null(ordering)) {
    # Validate ordering argument
    ordering <- match.arg(ordering, c("median", "mean", "lfdr"))

    order_result <- fash_post_ordering(eb_output_selected, ordering = ordering)
    posterior_weights_matrix <- order_result$ordered_matrix
    ordered_indices <- order_result$ordered_indices
  } else {
    ordered_indices <- seq_len(nrow(posterior_weights_matrix))
  }

  # Extract PSD values
  psd_values <- as.numeric(colnames(posterior_weights_matrix))

  # Convert the posterior matrix to a data frame for ggplot
  posterior_weights_df <- as.data.frame(posterior_weights_matrix)
  posterior_weights_df$id <- ordered_indices

  # Melt the data frame for ggplot
  melted_data <- reshape2::melt(posterior_weights_df, id.vars = "id")
  melted_data$variable <- as.numeric(as.character(melted_data$variable))

  # Adjust the PSD variable for discrete or continuous plotting
  if (discrete) {
    # Round PSD values and convert to factor
    melted_data$variable <- factor(round(melted_data$variable, 3), levels = round(psd_values, 3))
    fill_scale <- ggplot2::scale_fill_brewer(palette = "Set3", name = "PSD (Rounded)")
  } else {
    fill_scale <- ggplot2::scale_fill_gradient(low = "white", high = "blue", name = "PSD")
  }

  # Create the structure plot
  melted_data$id <- factor(melted_data$id,levels = posterior_weights_df$id)
  return(ggplot2::ggplot(melted_data,
                         ggplot2::aes(x = .data$id, y = .data$value,
                                      fill = .data$variable)) +
           ggplot2::geom_bar(stat = "identity", position = "stack") +
           ggplot2::labs(
             x = "Datasets",
             y = "Posterior Weight",
             title = "Structure Plot of Posterior Weights"
           ) +
           fill_scale +
           ggplot2::coord_flip() +
           ggplot2::theme_minimal() +
           ggplot2::theme(
             axis.text.y = ggplot2::element_blank(),
             axis.ticks.y = ggplot2::element_blank(),
             panel.grid = ggplot2::element_blank(),
             panel.background = ggplot2::element_rect(fill = "white"),
             plot.background = ggplot2::element_rect(fill = "white")
           ))
}



#' Heatmap Plot of Posterior Weights for FASH Objects
#'
#' This function generates a heatmap plot visualizing the posterior weights from a \code{fash} object.
#' The y-axis shows dataset names, the x-axis shows PSD grid values, and point sizes
#' represent the posterior weights.
#'
#' @param object A \code{fash} object containing posterior weights.
#' @param selected_indices Optional character vector of dataset names or numeric indices
#'   to specify which rows (datasets) to display. Default is \code{NULL} (all datasets).
#' @param size_range A numeric vector of length 2 specifying the range of point sizes. Default is \code{c(1, 8)}.
#' @param size_breaks A numeric vector specifying size breaks from 0.1 to 0.9.
#'   Default is \code{NULL}, which automatically selects a set of breaks.
#' @param font_size A numeric value specifying the base font size for theme elements. Default is \code{10}.
#' @param ... Additional arguments passed to \code{ggplot2::theme} or \code{ggplot2::geom_point}.
#'
#' @return A \code{ggplot} object representing the heatmap plot of posterior weights.
#'
#' @examples
#' # Simulate example
#' data_list <- lapply(1:10, function(i) data.frame(y = rpois(16, 5), x = 1:16, offset = 0))
#' grid <- seq(0, 2, length.out = 6)
#' fash_obj <- fash(data_list = data_list, Y = "y", smooth_var = "x", grid = grid, likelihood = "poisson")
#'
#' # Heatmap plot for all datasets
#' plot_heatmap(fash_obj)
#'
#' # Subset some datasets
#' plot_heatmap(fash_obj, selected_indices = 1:5)
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_size theme element_text
#' @importFrom cowplot theme_cowplot
#' @importFrom rlang .data
#'
#' @export
plot_heatmap <- function(object,
                         selected_indices = NULL,
                         size_range = c(1, 8),
                         size_breaks = NULL,
                         font_size = 10,
                         ...) {
  if (!inherits(object, "fash")) stop("Input must be a fash object.")

  posterior_weights <- object$posterior_weights

  # Optionally subset rows
  if (!is.null(selected_indices)) {
    if (is.character(selected_indices)) {
      posterior_weights <- posterior_weights[selected_indices, , drop = FALSE]
    } else if (is.numeric(selected_indices)) {
      posterior_weights <- posterior_weights[selected_indices, , drop = FALSE]
    } else {
      stop("selected_indices must be a character or numeric vector.")
    }
  }

  # Auto-generate size_breaks if NULL
  if (is.null(size_breaks)) {
    size_breaks <- c(0.1,0.3,0.5,0.7,0.9)
    size_breaks <- unique(size_breaks)  # Ensure unique breaks
  } else {
    size_breaks <- sort(unique(size_breaks))  # Ensure unique and sorted breaks
  }

  # Round the colname
  colnames(posterior_weights) <- round(as.numeric(colnames(posterior_weights)), 3)

  # Reshape for ggplot
  pdat <- data.frame(
    "dataset" = rep(rownames(posterior_weights),
                    times = ncol(posterior_weights)),
    "psd"     = rep(colnames(posterior_weights),
                    each = nrow(posterior_weights)),
    "weight"  = as.vector(posterior_weights)
  )
  pdat$dataset <- factor(pdat$dataset, levels = rev(unique(pdat$dataset)))
  pdat$psd     <- factor(pdat$psd, levels = unique(pdat$psd))

  # Make plot
  p <- ggplot2::ggplot(pdat, ggplot2::aes(x = .data$psd,
                                          y = .data$dataset,
                                          size = .data$weight)) +
    ggplot2::geom_point(shape = 21, fill = "black", color = "white", ...) +
    ggplot2::scale_size(range = size_range, breaks = size_breaks) +
    cowplot::theme_cowplot(font_size = font_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

  return(p)
}


#' Simulate functions from an integrated Wiener process prior
#'
#' This internal helper simulates sample paths from an integrated Wiener process
#' (IWP)–type prior using a spline basis plus a low-order global polynomial
#' trend. The spline coefficients are drawn from a multivariate normal
#' distribution with a precision matrix constructed from the knot locations,
#' scaled to achieve a target predictive standard deviation (PSD). The global
#' polynomial coefficients are drawn independently from a normal prior.
#'
#' @param n_samps Integer; number of function samples to draw.
#' @param n_basis Integer; number of spline basis functions (i.e., number of
#'   knots). Must be at least 3.
#' @param psd Numeric; target predictive standard deviation (PSD) for the IWP
#'   component, on the scale implied by \code{pred_step}.
#' @param sd_poly Numeric; standard deviation for the global polynomial
#'   coefficients.
#' @param p Integer; order of the integrated Wiener process / local polynomial.
#' @param pred_step Numeric; prediction step size used in the PSD scaling
#'   formula.
#' @param x_range Numeric vector of length 2 giving the range of the input
#'   domain \code{[x_min, x_max]}. If \code{NULL}, defaults to \code{c(0, 10)}.
#'   Used to construct the knot locations.
#' @param x_new Optional numeric vector giving the evaluation points at which
#'   the sampled functions should be returned. If \code{NULL}, 100 equally
#'   spaced points over \code{x_range} are used.
#'
#' @details
#' The function constructs a spline basis over \code{x_range} using
#' \code{n_basis} knots, and a global polynomial design of order \code{p} over
#' \code{x_new}. Spline coefficients are sampled from a multivariate normal
#' distribution with precision matrix
#' \deqn{Q = \frac{1}{\sigma^2} Q_{\text{IWP}}(knots),}
#' where \code{Q_IWP} is computed by \code{compute_weights_precision_helper()}
#' and \code{sigma} is chosen so that the resulting process has predictive
#' standard deviation approximately equal to \code{psd} at step size
#' \code{pred_step}. The global polynomial coefficients are sampled
#' independently from \eqn{N(0, \text{sd\_poly}^2)}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{samples}{A numeric matrix of dimension
#'     \code{length(x_new) × n_samps} containing the simulated function values.
#'     Each column corresponds to one function sample, and each row corresponds
#'     to an evaluation point in \code{x_new}.}
#'   \item{x_new}{The numeric vector of evaluation points at which the
#'     functions are simulated.}
#' }
#'
#' @examples
#' \dontrun{
#'   sim <- simulate_IWP(n_samps = 5, n_basis = 20, psd = 1, p = 2)
#'   # Each column is a sample
#'   matplot(sim$x_new, sim$samples, type = "l", lty = 1,
#'           xlab = "x", ylab = "Sampled function")
#' }
#'
#' @keywords internal
simulate_IWP <- function(n_samps = 1,
                         n_basis = 20,
                         psd = 1,
                         sd_poly = 0.1,
                         p = 1,
                         pred_step = 1,
                         x_range = NULL,
                         x_new = NULL) {
  if (n_basis < 3L) {
    stop("n_basis must be greater than 3.")
  }
  if (p < 1L) {
    stop("p must be a positive integer.")
  }

  if (is.null(x_range)) {
    x_range <- c(0, 10)
  }

  if (is.null(x_new)) {
    x_new <- seq(x_range[1], x_range[2], length.out = 100)
  }

  x_min <- x_range[1]
  x_max <- x_range[2]

  # Optional: warn if extrapolating outside the knot range
  if (any(x_new < x_min | x_new > x_max)) {
    warning("Some x_new values lie outside x_range; extrapolation may be unreliable.")
  }

  # Knots for the spline basis
  knots <- seq(x_min, x_max, length.out = n_basis)

  spline_new   <- fashr:::local_poly_helper(knots = knots, refined_x = x_new, p = p)
  x_new_design <- fashr:::global_poly_helper(x = x_new, p = p)

  if(psd > 0){
    # Precision for spline weights (IWP part), scaled to match target PSD
    sd_function <- psd / sqrt(
      (pred_step^((2 * p) - 1)) /
        (((2 * p) - 1) * (factorial(p - 1)^2))
    )

    prec_mat <- (1 / sd_function^2) * fashr:::compute_weights_precision_helper(knots)

    # Basis weights for the IWP component
    weights <- LaplacesDemon::rmvnp(
      n    = n_samps,
      mu   = rep(0, ncol(prec_mat)),
      Omega = prec_mat
    )
  }
  else{
    weights <- matrix(0, nrow = n_samps, ncol = ncol(spline_new))
  }

  # Global polynomial component
  beta <- rnorm(n = p * n_samps, mean = 0, sd = sd_poly)
  beta_matrix <- matrix(beta, nrow = n_samps, ncol = p, byrow = TRUE)


  # Combine global polynomial and IWP components
  samps <- beta_matrix %*% t(x_new_design) + weights %*% t(spline_new)
  # Transpose so that columns are samples and rows are x_new
  samps <- t(samps)

  list(
    samples = samps,
    x_new   = x_new
  )
}



#' Simulate sample paths from the FASH prior
#'
#' This function simulates sample paths from the prior over effect functions
#' implied by a fitted \code{fash} object. The prior is treated as a finite
#' mixture over predictive standard deviation (PSD) values stored in
#' \code{fash_obj$prior_weights}. Each sample path is drawn from an integrated
#' Wiener process (IWP) prior plus a global polynomial trend.
#'
#' The function provides optional constraints on the global polynomial
#' component:
#' \itemize{
#'   \item \code{"none"}: use the polynomial variance from the fitted model (default)
#'   \item \code{"initial"}: force the polynomial part to be identically zero
#'         (equivalent to letting betaprec = Inf)
#'   \item \code{"orthogonal"}: regress out the polynomial trend from each sample
#'         so that the resulting sample path is orthogonal to all global
#'         polynomial basis functions
#' }
#'
#' @param fash_obj A fitted \code{fash} object containing:
#'   \itemize{
#'     \item \code{prior_weights}: a data frame with columns \code{psd} and \code{prior_weight};
#'     \item \code{settings}: a list with elements \code{num_basis}, \code{order},
#'           \code{pred_step}, \code{betaprec}.
#'   }
#' @param M Integer; total number of prior samples to draw.
#' @param constraints Character; one of:
#'   \code{"none"}, \code{"initial"}, \code{"orthogonal"}.
#' @param x_range Optional numeric vector of length 2 defining the simulation
#'   domain. If missing, inferred from the data.
#' @param x_new Optional numeric vector giving evaluation points for the samples.
#'
#' @return A list containing:
#' \describe{
#'   \item{samples}{A matrix of size \code{length(x_new) × M}; each column is a sampled function.}
#'   \item{x_new}{Evaluation grid.}
#'   \item{psd}{Length-\code{M} vector giving the PSD used for each sample.}
#'   \item{component}{Length-\code{M} vector giving mixture component index used.}
#'   \item{prior_weights}{Original prior weights.}
#'   \item{settings}{Settings used in the simulation.}
#' }
#'
#' @export
simulate_fash_prior <- function(fash_obj,
                                M = 100,
                                constraints = c("none", "initial", "orthogonal"),
                                x_range = NULL,
                                x_new = NULL) {

  constraints <- match.arg(constraints)

  # ------------------------------------------------------------
  # 1. Extract prior weights
  # ------------------------------------------------------------
  prior_df <- fash_obj$prior_weights
  if (is.null(prior_df)) stop("`fash_obj$prior_weights` is NULL.")
  if (!all(c("psd", "prior_weight") %in% names(prior_df)))
    stop("`prior_weights` must contain columns `psd` and `prior_weight`.")

  # normalize weights
  w <- prior_df$prior_weight
  w <- w / sum(w)
  K <- nrow(prior_df)

  # ------------------------------------------------------------
  # 2. Extract settings
  # ------------------------------------------------------------
  settings <- fash_obj$settings
  if (is.null(settings))
    stop("`fash_obj$settings` is NULL.")

  n_basis   <- settings$num_basis
  p         <- settings$order
  pred_step <- settings$pred_step
  betaprec  <- settings$betaprec

  # constraints = "initial" → polynomial ≡ 0
  if (constraints == "initial") {
    betaprec <- Inf
  }

  # compute sd_poly
  sd_poly <- if (betaprec > 0) 1 / sqrt(betaprec) else 0

  # ------------------------------------------------------------
  # 3. Infer x_range if needed
  # ------------------------------------------------------------
  if (is.null(x_range) && is.null(x_new)) {
    data_list <- fash_obj$fash_data$data_list
    if (!is.null(data_list) && length(data_list) > 0L)
      x_range <- range(data_list[[1]]$x)
  }

  # ------------------------------------------------------------
  # 4. Mixture sampling
  # ------------------------------------------------------------
  comp_idx <- sample.int(K, size = M, replace = TRUE, prob = w)
  uniq_comp <- sort(unique(comp_idx))

  x_new_common <- NULL
  samples_mat  <- NULL

  # ------------------------------------------------------------
  # 5. Generate samples component-by-component
  # ------------------------------------------------------------
  for (u in uniq_comp) {
    n_u <- sum(comp_idx == u)

    sim_u <- simulate_IWP(
      n_samps   = n_u,
      n_basis   = n_basis,
      psd       = prior_df$psd[u],
      sd_poly   = sd_poly,
      p         = p,
      pred_step = pred_step,
      x_range   = x_range,
      x_new     = x_new
    )

    if (is.null(x_new_common)) {
      x_new_common <- sim_u$x_new
      samples_mat  <- matrix(NA_real_,
                             nrow = length(x_new_common),
                             ncol = M)
    }

    samples_mat[, comp_idx == u] <- sim_u$samples
  }

  # ------------------------------------------------------------
  # 6. constraints = "orthogonal":
  #    regress out polynomial component from each sample
  # ------------------------------------------------------------
  if (constraints == "orthogonal") {
    X_poly <- fashr:::global_poly_helper(x = x_new_common, p = p)  # n × p
    XtX <- crossprod(X_poly)
    XtX_inv <- solve(XtX)
    Xt <- t(X_poly)

    for (j in seq_len(M)) {
      f_j <- samples_mat[, j]
      beta_hat <- XtX_inv %*% (Xt %*% f_j)
      fitted_poly <- X_poly %*% beta_hat
      samples_mat[, j] <- f_j - fitted_poly
    }
  }

  list(
    samples       = samples_mat,
    x_new         = x_new_common,
    psd           = prior_df$psd[comp_idx],
    component     = comp_idx,
    prior_weights = prior_df,
    settings      = settings
  )
}






#' Visualize the FASH prior over effect functions
#'
#' This function visualizes the fitted fash prior in the
#' \code{fash} object, using simulation via \code{simulate_fash_prior()}.
#' Two types of visualizations are supported:
#' \itemize{
#'   \item \code{"sample_path"}: plot individual sample paths from the prior,
#'         with different PSD mixture components shown in different colors and
#'         line types.
#'   \item \code{"psd"}: plot a histogram of the PSD values corresponding to
#'         the \code{M} prior samples drawn from the mixture.
#' }
#'
#' @param fash_obj A fitted \code{fash} object, passed to
#'   \code{simulate_fash_prior()}.
#' @param plot_type Character string specifying the type of plot:
#'   \code{"sample_path"} or \code{"psd"}.
#' @param M Integer; number of prior samples to draw via
#'   \code{simulate_fash_prior()}. For \code{plot_type = "sample_path"}, this is
#'   the number of sample paths drawn (default \code{M = 100}). For
#'   \code{plot_type = "psd"}, this is the number of PSD values shown in the
#'   histogram (i.e., \code{M} PSD values simulated from the prior mixture).
#' @param constraints Character; passed to \code{simulate_fash_prior()}.
#'   One of \code{"none"}, \code{"initial"}, or \code{"orthogonal"}.
#' @param x_range Optional numeric vector of length 2 defining the simulation
#'   domain. If \code{NULL}, may be inferred inside \code{simulate_fash_prior()}.
#' @param x_new Optional numeric vector of evaluation points for the prior
#'   functions. If \code{NULL}, a default grid is used inside
#'   \code{simulate_fash_prior()}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
visualize_fash_prior <- function(fash_obj,
                                 plot_type = c("sample_path", "psd"),
                                 M = 100,
                                 constraints = c("none", "initial", "orthogonal"),
                                 x_range = NULL,
                                 x_new = NULL,
                                 ...) {

  plot_type   <- match.arg(plot_type)
  constraints <- match.arg(constraints)

  # 1. Simulate from prior
  sim <- simulate_fash_prior(
    fash_obj    = fash_obj,
    M           = M,
    constraints = constraints,
    x_range     = x_range,
    x_new       = x_new
  )

  x_grid   <- sim$x_new
  samples  <- sim$samples  # matrix: length(x_grid) x M
  psd_vec  <- sim$psd      # length M, PSD per sample

  if (plot_type == "sample_path") {

    df <- data.frame(
      x         = rep(x_grid, times = M),
      y         = as.vector(samples),
      sample_id = rep(seq_len(M), each = length(x_grid)),
      psd       = factor(round(rep(psd_vec, each = length(x_grid)), 4))
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y,
                                          group = sample_id,
                                          color = psd,
                                          linetype = psd)) +
      ggplot2::geom_line(alpha = 0.6) +
      ggplot2::labs(
        x = "Condition",
        y = "Function value",
        color    = "PSD",
        linetype = "PSD",
        title = "Sample paths from FASH prior"
      ) +
      ggplot2::theme_minimal()

    # symmetric y-limits based on overall magnitude
    max_abs <- max(abs(samples), na.rm = TRUE)
    if (is.finite(max_abs) && max_abs > 0) {
      radius <- max_abs * 1.05
      p <- p + ggplot2::coord_cartesian(ylim = c(-radius, radius))
    }

    return(p)

  } else if (plot_type == "psd") {

    df_psd <- data.frame(psd = psd_vec)

    p <- ggplot2::ggplot(df_psd, ggplot2::aes(x = psd)) +
      ggplot2::geom_histogram(
        bins = min(30, max(5, length(unique(psd_vec)))),
        color = "white"
      ) +
      ggplot2::labs(
        x = "Predictive standard deviation (PSD)",
        y = "Frequency",
        title = sprintf("Histogram of %d PSD values from FASH prior", M)
      ) +
      ggplot2::theme_minimal()

    return(p)
  }
}







