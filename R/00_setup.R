#' Roxygen commands for package registration
#'
#' @useDynLib Gaussian_dep
#' @useDynLib Gaussian_ind
#' @useDynLib Gaussian_dep_fixed
#' @useDynLib Gaussian_ind_fixed
#' @useDynLib Poisson_ind
#' @useDynLib Poisson_ind_fixed
#' @importFrom TMB compile
#'
#' @keywords internal
#'
dummy <- function() {
  NULL
}

#' Validate the requested number of cores for parallel processing
#'
#' \code{parallel::mclapply} relies on forking, which is not available on
#' Windows. This helper falls back to sequential execution (with a warning)
#' on platforms where forking is unsupported.
#'
#' @param num_cores An integer specifying the requested number of cores.
#' @return An integer number of cores that is safe to use on this platform.
#'
#' @keywords internal
#'
sanitize_num_cores <- function(num_cores) {
  if (!is.numeric(num_cores) || length(num_cores) != 1 || num_cores < 1 || num_cores %% 1 != 0) {
    stop("num_cores must be a positive integer.")
  }
  if (num_cores > 1 && .Platform$OS.type == "windows") {
    warning("Parallel processing via forking is not supported on Windows; using num_cores = 1.")
    return(1L)
  }
  as.integer(num_cores)
}

#' Constructing and evaluating the global polynomials, to account for boundary conditions (design matrix)
#'
#' @param x A vector of locations to evaluate the global polynomials
#' @param p An integer value indicates the order of smoothness
#' @return A matrix with i,j component being the value of jth basis function
#' value at ith element of x, the ncol should equal to p, and nrow
#' should equal to the number of elements in x
#' @examples
#' fashr:::global_poly_helper(x = c(0, 0.2, 0.4, 0.6, 0.8), p = 2)
#'
#' @keywords internal
#'
global_poly_helper <- function(x, p = 2) {
  result <- NULL
  for (i in 1:p) {
    result <- cbind(result, x^(i - 1))
  }
  result
}

get_local_poly <- function(knots, refined_x, p) {
  dif <- diff(knots)
  nn <- length(refined_x)
  n <- length(knots)
  D <- matrix(0, nrow = nn, ncol = n - 1)
  for (j in 1:nn) {
    for (i in 1:(n - 1)) {
      if (refined_x[j] <= knots[i]) {
        D[j, i] <- 0
      } else if (refined_x[j] <= knots[i + 1] & refined_x[j] >= knots[i]) {
        D[j, i] <- (1 / factorial(p)) * (refined_x[j] - knots[i])^p
      } else {
        k <- 1:p
        D[j, i] <- sum((dif[i]^k) * ((refined_x[j] - knots[i + 1])^(p - k)) / (factorial(k) * factorial(p - k)))
      }
    }
  }
  D # Local poly design matrix
}

#' Constructing and evaluating the local O-spline basis (design matrix)
#'
#' @param knots A vector of knots used to construct the O-spline basis, first knot should be viewed as "0",
#' the reference starting location. These k knots will define (k-1) basis function in total.
#' @param refined_x A vector of locations to evaluate the O-spline basis
#' @param p An integer value indicates the order of smoothness
#' @param neg_sign_order An integer value N such that D = ((-1)^N)*D for the splines at negative knots.
#' @return A matrix with i,j component being the value of jth basis function
#' value at ith element of refined_x, the ncol should equal to number of knots minus 1, and nrow
#' should equal to the number of elements in refined_x.
#' @examples
#' fashr:::local_poly_helper(knots = c(0, 0.2, 0.4, 0.6, 0.8), refined_x = seq(0, 0.8, by = 0.1), p = 2)
#'
#' @keywords internal
#'
local_poly_helper <- function(knots, refined_x, p = 2, neg_sign_order = 0) {
  if (min(knots) >= 0) {
    # The case of all-positive knots
    D <- get_local_poly(knots, refined_x, p)
  } else if (max(knots) <= 0) {
    # Handle the negative part only
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    D <- get_local_poly(knots_neg, refined_x_neg, p)
    D <- D * ((-1)^neg_sign_order)
  } else {
    # Handle the negative part
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    D1 <- get_local_poly(knots_neg, refined_x_neg, p)
    D1 <- D1 * ((-1)^neg_sign_order)

    # Handle the positive part
    refined_x_pos <- ifelse(refined_x > 0, refined_x, 0)
    knots_pos <- unique(sort(ifelse(knots > 0, knots, 0)))
    D2 <- get_local_poly(knots_pos, refined_x_pos, p)

    D <- cbind(D1, D2)
  }
  D # Local poly design matrix
}

#' Constructing the precision matrix given the knot sequence (helper)
#'
#' @param x A vector of knots used to construct the O-spline basis, first knot should be viewed as "0",
#' the reference starting location. These k knots will define (k-1) basis function in total.
#' @return A precision matrix of the corresponding basis function, should be diagonal matrix with
#' size (k-1) by (k-1).
#' @examples
#' fashr:::compute_weights_precision_helper(x = c(0,0.2,0.4,0.6,0.8))
#'
#' @keywords internal
#'
compute_weights_precision_helper <- function(x){
  d <- diff(x)
  Precweights <- diag(d)
  Precweights
}

#' Setup Quantities for FASH
#'
#' This function prepares the necessary data structures for implementing the Function Adaptive SHrinkage (FASH) approach. It processes the input data, smoothing variables, offsets, standard errors, and optional precision matrices to generate a list of data frames and supporting matrices.
#'
#' @param Y Either a numeric matrix of response variables or a character string specifying the column name in \code{data_list} for response variables. Required if \code{data_list} is provided.
#' @param smooth_var A numeric matrix, vector, or a character string specifying the column name in \code{data_list} for smoothing variables. Required if \code{data_list} is provided.
#' @param offset A numeric matrix, vector, scalar, or a character string specifying the column name in \code{data_list} for offset variables.
#' @param S A numeric matrix, vector, scalar, or list representing the standard errors of \code{Y}. Or a character string specifying the column name in \code{data_list} for SD. If provided as a matrix, it must have the same dimensions as \code{Y}. If a vector or scalar, it will be expanded to match the dimensions of \code{Y}.
#' @param Omega Either a list of precision matrices (one for each dataset) or a single precision matrix (shared across all datasets). Default is \code{NULL}.
#' @param data_list A list of data frames, where each data frame corresponds to a single dataset.
#'
#' @return A list with the following components:
#' \item{data_list}{A list of data frames, where each data frame corresponds to a row of \code{Y}. Each data frame contains:
#'   \describe{
#'     \item{\code{y}}{The response variables for the corresponding row of \code{Y}.}
#'     \item{\code{x}}{The smoothing variables for the corresponding row of \code{smooth_var}.}
#'     \item{\code{offset}}{The offset values for the corresponding row of \code{offset}.}
#'   }
#' }
#' \item{S}{A list of standard errors, where each element corresponds to the standard errors for a single dataset.}
#' \item{Omega}{A list of precision matrices, where each element corresponds to a single dataset's precision matrix.}
#'
#' @examples
#' # Example usage with matrix input
#' set.seed(1)
#' Y <- matrix(rnorm(20), nrow = 4, ncol = 5)
#' smooth_var <- matrix(runif(20), nrow = 4, ncol = 5)
#' offset <- 1
#' S <- c(0.5, 0.8, 1.2, 1.0, 0.9)
#' Omega <- diag(c(1, 2, 3, 4, 5))
#' result <- fash_set_data(Y = Y, smooth_var = smooth_var, offset = offset, S = S, Omega = Omega)
#'
#' # Example usage with data_list input
#' data_list <- list(
#'   data.frame(y = rnorm(5), x = 1:5, offset = 0),
#'   data.frame(y = rnorm(5), x = 1:5, offset = 0)
#' )
#' S_list <- list(rep(0.5, 5), rep(0.8, 5))
#' Omega_list <- list(diag(5), diag(5))
#' result <- fash_set_data(data_list = data_list, Y = "y", smooth_var = "x", offset = "offset", S = S_list, Omega = Omega_list)
#'
#' @export
#' 
fash_set_data <- function(Y, smooth_var, offset = 0, S = NULL, Omega = NULL, data_list = NULL) {
  # If data_list is provided, extract variables
  if (!is.null(data_list)) {
    if (!is.list(data_list) || !all(sapply(data_list, is.data.frame))) {
      stop("data_list must be a list of data frames.")
    }
    n_datasets <- length(data_list)

    # Extract variables from data_list
    processed_data_list <- lapply(data_list, function(df) {
      # Extract response variable (Y)
      if (!Y %in% names(df)) stop(sprintf("Variable '%s' not found in data frame.", Y))
      y <- df[[Y]]

      # Extract smoothing variable (smooth_var)
      if (!smooth_var %in% names(df)) stop(sprintf("Variable '%s' not found in data frame.", smooth_var))
      x <- df[[smooth_var]]

      # Extract or process offset
      if (is.character(offset)) {
        if (!offset %in% names(df)) {
          stop(sprintf("Offset variable '%s' not found in data frame.", offset))
        }
        offset_value <- df[[offset]]
      } else if (is.numeric(offset) && length(offset) == 1) {
        offset_value <- rep(offset, length(y))
      } else {
        stop("When using data_list, offset must be a valid column name or a scalar value.")
      }

      # Return processed data frame
      data.frame(y = y, x = x, offset = offset_value)
    })

  } else {
    # If data_list is NULL, construct it from Y, smooth_var, and offset
    n_datasets <- nrow(Y)

    # Validate and process smooth_var
    if (is.matrix(smooth_var)) {
      if (nrow(smooth_var) != n_datasets || ncol(smooth_var) != ncol(Y)) {
        stop("smooth_var must have the same dimensions as Y if it is a matrix.")
      }
    } else if (is.vector(smooth_var) && length(smooth_var) == ncol(Y)) {
      smooth_var <- matrix(rep(smooth_var, n_datasets), nrow = n_datasets, byrow = TRUE)
    } else {
      stop("smooth_var must be a matrix or a vector matching the number of columns in Y.")
    }

    # Validate and process offset
    if (is.matrix(offset)) {
      if (nrow(offset) != n_datasets || ncol(offset) != ncol(Y)) {
        stop("offset must have the same dimensions as Y if it is a matrix.")
      }
    } else if (is.vector(offset) && length(offset) == ncol(Y)) {
      offset <- matrix(rep(offset, n_datasets), nrow = n_datasets, byrow = TRUE)
    } else if (length(offset) == 1) {
      offset <- matrix(offset, nrow = n_datasets, ncol = ncol(Y))
    } else {
      stop("offset must be a matrix, vector, or scalar.")
    }

    # Create processed_data_list
    processed_data_list <- lapply(1:n_datasets, function(i) {
      data.frame(y = Y[i,], x = smooth_var[i,], offset = offset[i,])
    })
  }

  # Validate and process S
  if (!is.null(S)) {
    if (is.list(S)) {
      if (length(S) != n_datasets) {
        stop("S must have the same length as the number of datasets when provided as a list.")
      }
    } else if (is.matrix(S)) {
      S <- split(S, row(S))
    } else if (is.vector(S) & is.numeric(S)) {
      S <- split(matrix(rep(S, n_datasets), nrow = n_datasets, byrow = TRUE), 1:n_datasets)
    } else if(is.character(S)){
      # check if S is a column name in data_list
      if(is.null(data_list)){
        stop("data_list must be provided if S is a character string.")
      }
      if(!all(S %in% names(data_list[[1]]))){
        stop("S must be a valid column name in data_list.")
      }
      S <- lapply(data_list, function(df) df[[S]])
    }
    else {
      stop("S must be a string, list, matrix, or numeric vector.")
    }
  }

  # Validate and process Omega
  if (!is.null(Omega)) {
    if (is.list(Omega)) {
      if (length(Omega) != n_datasets) {
        stop("Omega must have the same length as the number of datasets when provided as a list.")
      }
    } else if (is.matrix(Omega)) {
      Omega <- replicate(n_datasets, Omega, simplify = FALSE)
    } else {
      stop("Omega must be either a list of matrices or a single matrix.")
    }
  }

  # Return structured list
  return(list(data_list = processed_data_list, S = S, Omega = Omega))
}

#' Create a TMB Data Object for a Specific Dataset
#'
#' Processes the \code{tmbdat} object for a specific dataset, using the provided smoothing variables,
#' standard errors, and precision matrix. This function generates spline basis matrices,
#' penalty matrices, and other quantities required for modeling.
#'
#' @param data_i A single dataset extracted from the \code{data_list} component of \code{fash_set_data}.
#'               Must be a list containing \code{y}, \code{x}, and \code{offset}.
#' @param Si A numeric vector representing the standard errors for the dataset.
#' @param Omegai A numeric precision matrix for the dataset.
#' @param num_basis An integer specifying the number of O-Spline basis functions to use for the approximation.
#' @param betaprec A numeric value representing the precision of the fixed effects coefficients.
#' @param order An integer specifying the order of the Integrated Wiener Process (IWP) prior.
#'
#' @return A list formatted as a \code{tmbdat} object, containing:
#' \item{y}{A numeric vector of response variables for the dataset.}
#' \item{X}{A sparse matrix representing the design matrix for fixed effects.}
#' \item{P}{A sparse matrix representing the penalty matrix for the spline coefficients.}
#' \item{B}{A sparse matrix representing the design matrix for random effects.}
#' \item{offset}{A numeric vector of offsets for the dataset.}
#' \item{logPdet}{The log determinant of the penalty matrix.}
#' \item{betaprec}{The precision for the fixed effects coefficients.}
#' \item{S}{The standard errors for the dataset, if provided as \code{Si}.}
#' \item{Omega}{The precision matrix for the dataset as a sparse matrix, if provided as \code{Omegai}.}
#' \item{log_det_Omega}{The log determinant of the precision matrix, if \code{Omegai} is provided.}
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
#' tmbdat <- fash_set_tmbdat(data$data_list[[1]], Si = S[[1]], Omegai = Omega)
#'
#' @importFrom methods as
#'
#' @export
#'
fash_set_tmbdat <- function(data_i, Si = NULL, Omegai = NULL, num_basis = 30, betaprec = 1e-6, order = 1) {
  # Extract smoothing variables and response
  y <- data_i$y
  x <- data_i$x
  offset <- data_i$offset

  # Generate spline knots for the smoothing variable
  knots <- seq(min(x), max(x), length.out = num_basis)

  # Generate the design matrices for the spline basis
  X <- global_poly_helper(x, p = order)
  P <- compute_weights_precision_helper(knots)
  B <- local_poly_helper(knots = knots, refined_x = x, p = order)

  # Prepare the tmbdat list
  tmbdat <- list(
    y = y,                                                    # Response variable
    X = as(as(as(X, "dMatrix"), "generalMatrix"), "TsparseMatrix"), # Design matrix for fixed effects
    P = as(as(as(P, "dMatrix"), "generalMatrix"), "TsparseMatrix"), # Penalty matrix for spline coefficients
    B = as(as(as(B, "dMatrix"), "generalMatrix"), "TsparseMatrix"), # Design matrix for random effects
    offset = offset,                                          # Offset term
    logPdet = as.numeric(determinant(P, logarithm = TRUE)$modulus), # Log determinant of penalty matrix
    betaprec = betaprec                                       # Precision for fixed effects
  )

  # Add Omega and compute log determinant if available
  if (!is.null(Si)) {
    # check if Si has the same same length as y
    if (length(Si) != length(y)) {
      # if not, see if Si is a scalar
      if (length(Si) == 1) {
        Si <- rep(Si, length(y))
      } else{
        stop("Standard errors (Si) must have the same length as the response variable (y).")
      }
    }
    tmbdat$S <- Si
  }
  if (!is.null(Omegai)) {
    # check if the dimension of Omegai matches the number of observations
    if (nrow(Omegai) != length(y) || ncol(Omegai) != length(y)) {
      stop("Precision matrix (Omegai) must have the same dimensions as the response variable (y).")
    }
    tmbdat$Omega <- as(as(as(Omegai, "dMatrix"), "generalMatrix"), "TsparseMatrix")
    tmbdat$log_det_Omega <- as.numeric(determinant(Omegai, logarithm = TRUE)$modulus)
  }

  return(tmbdat)
}

