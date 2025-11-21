#' Simulate Observed Data from a Specified True Function
#'
#' This function generates synthetic observed data by evaluating a user-specified
#' true function \code{g} on a grid of x-values and adding Gaussian noise
#' with heteroscedastic standard deviation sampled from a provided vector or scalar.
#'
#' @param g A function specifying the true underlying signal to simulate.
#' @param x A numeric vector of grid points where the function \code{g} is evaluated.
#' @param sd A numeric value or vector specifying the standard deviation(s) of noise.
#'   If a scalar, all points have the same standard deviation. If a vector, it is
#'   sampled with replacement for each \code{x}.
#'
#' @return A data frame containing:
#'   \describe{
#'     \item{\code{x}}{Input grid points.}
#'     \item{\code{y}}{Simulated noisy observations.}
#'     \item{\code{truef}}{True function values at \code{x}.}
#'     \item{\code{sd}}{Standard deviation used for noise at each point.}
#'   }
#'
#' @examples
#' g <- function(x) sin(x / 5)
#' fashr:::simulate_data(g, sd = 0.2)
#'
#' @importFrom stats rnorm
#'
#' @keywords internal
#'
simulate_data <- function(g, x = NULL, sd = 0.1){
  if (is.null(x)) {
    x <- 1:16
  }
  # simulate sd from sampling from sd with replacement
  sd <- sample(x = sd, size = length(x), replace = TRUE)
  y <- g(x) + rnorm(n = length(x), sd = sd, mean = 0)
  return(data.frame(x = x, y = y, truef = g(x), sd = sd))
}

#' Simulate a Random Nonlinear Function from an Integrated Wiener Process (IWP) Prior
#'
#' This function generates a random smooth function sampled from an Integrated Wiener Process (IWP) prior,
#' parameterized by the number of basis functions and smoothing parameters.
#'
#' @param n_basis An integer specifying the number of spline basis functions (minimum 3).
#' @param sd_function A numeric value representing the predictive standard deviation of the function.
#' @param sd_poly A numeric value specifying the standard deviation for the linear polynomial component.
#' @param p An integer specifying the order of the polynomial trend component (default is 1 for linear).
#' @param pred_step A numeric value defining the prediction step size that determines the predictive standard deviation.
#' @param x_range A numeric vector of length 2 specifying the range of x-values over which the function is defined.
#'
#' @return A function that can evaluate the generated random smooth signal at any new \code{x} values:
#'   \describe{
#'     \item{function}{Returns the evaluated values of the sampled nonlinear function at specified \code{x} inputs.}
#'   }
#'
#' @examples
#' f <- fashr:::simulate_nonlinear_function()
#' plot(1:16, f(1:16), type = 'l')
#'
#' @importFrom stats rnorm
#' @importFrom LaplacesDemon rmvnp
#'
#' @keywords internal
#'
simulate_nonlinear_function <- function(n_basis = 20, sd_function = 1, sd_poly = 0.1, p = 1, pred_step = 16, x_range = NULL) {
  if(n_basis < 3) stop("n_basis must be greater than 3")

  if(is.null(x_range)) {
    x_range <- c(0, 16)
  }

  # Define the range and knots for the spline basis
  x_min <- x_range[1]
  x_max <- x_range[2]

  # Generate equally spaced knots within the range [x_min, x_max]
  knots <- seq(x_min, x_max, length.out = n_basis - 3)

  # Generate random weights for the basis functions
  sd_function <- sd_function/sqrt((pred_step^((2 * p) - 1)) / (((2 * p) - 1) * (factorial(p - 1)^2)))
  prec_mat <- (1/sd_function^2) * fashr:::compute_weights_precision_helper(knots)
  weights <- as.vector(LaplacesDemon::rmvnp(n = 1, mu = rep(0, ncol(prec_mat)), Omega = prec_mat))
  # Generate random weights for the linear functions
  beta_vec <- rnorm(n = p, mean = 0, sd = sd_poly)

  # Return a function that evaluates the spline at new x values
  function(x_new) {
    # Create the spline basis for the new x values using the predefined knots
    spline_new <- fashr:::local_poly_helper(knots = knots, refined_x = x_new, p = p)
    x_new_design <- fashr:::global_poly_helper(x = x_new, p = p)
    # Return the function
    return(x_new_design %*% beta_vec + as.vector(spline_new %*% weights))
  }
}



#' Simulate a Random Linear Function
#'
#' Generates a random linear function with intercept and slope drawn from normal distributions.
#'
#' @param sd_poly Standard deviation for both intercept and slope terms.
#' @param pred_step Resolution scaling for the slope variance.
#'
#' @return A function that evaluates the random linear function at new x-values.
#'
#' @examples
#'
#' f <- fashr:::simulate_linear_function()
#' plot(1:16, f(1:16), type = 'l')
#'
#' @importFrom stats rnorm
#'
#' @keywords internal
#'
simulate_linear_function <- function(sd_poly = 1, pred_step = 16){
  beta0 <- rnorm(1, mean = 0, sd = sd_poly)
  beta1 <- rnorm(1, mean = 0, sd = sd_poly/pred_step)
  function(x_new) {
    return(beta0 + beta1 * x_new)
  }
}

#' Simulate a Random Quadratic Function
#'
#' Generates a random quadratic function with coefficients drawn from normal distributions.
#'
#' @param sd_poly Standard deviation for the coefficients.
#'
#' @return A function that evaluates the random quadratic function at new x-values.
#'
#' @examples
#'
#' f <- fashr:::simulate_quadratic_function()
#' plot(1:16, f(1:16), type = 'l')
#'
#' @importFrom stats rnorm
#'
#' @keywords internal
#'
simulate_quadratic_function <- function(sd_poly = 1){
  beta0 <- rnorm(1, mean = 0, sd = sd_poly)
  beta1 <- rnorm(1, mean = 0, sd = sd_poly)
  beta2 <- rnorm(1, mean = 0, sd = sd_poly)
  function(x_new) {
    return(beta0 + beta1 * x_new + beta2 * x_new^2)
  }
}

#' Simulate a Random Constant Function (Nondynamic)
#'
#' Generates a constant function with value drawn from a normal distribution.
#'
#' @param sd_poly Standard deviation for the constant term.
#'
#' @return A function that evaluates to a constant at any x.
#'
#' @examples
#'
#' f <- fashr:::simulate_nondynamic_function()
#' plot(1:16, f(1:16), type = 'l')
#'
#' @importFrom stats rnorm
#'
#' @keywords internal
#'
simulate_nondynamic_function <- function(sd_poly = 1){
  beta0 <- rnorm(1, mean = 0, sd = sd_poly)
  Vectorize(function(x_new) {
    return(beta0)
  })
}

#' Simulate an Entire Observed Dataset from a Randomly Generated Function
#'
#' This function automates the process of generating a random underlying function
#' (of various types) and then sampling noisy observations from it, returning a complete simulated dataset.
#'
#' @param x A numeric vector specifying the grid points where the function is evaluated.
#' @param n_basis An integer specifying the number of basis functions (for \code{"nonlinear"} type).
#' @param sd_fun A numeric value controlling the smoothness of the function under the IWP prior.
#' @param sd A numeric value or vector specifying the standard deviation(s) of observation noise.
#' @param sd_poly A numeric value specifying the standard deviation for polynomial coefficients.
#' @param type A character string specifying the type of function to simulate. One of \code{"linear"}, \code{"quadratic"}, \code{"nonlinear"}, or \code{"nondynamic"}.
#' @param p An integer specifying the order of the polynomial trend in the nonlinear model.
#' @param pred_step A numeric value representing the resolution parameter for scaling the IWP precision.
#' @param normalize A logical value indicating whether to rescale the function to roughly zero mean and unit range.
#'
#' @return A data frame containing:
#'   \describe{
#'     \item{\code{x}}{Grid points where the function is evaluated.}
#'     \item{\code{y}}{Noisy observed values.}
#'     \item{\code{truef}}{True underlying function values at each \code{x}.}
#'     \item{\code{sd}}{Noise standard deviation used at each point.}
#'   }
#'
#' @examples
#' dat <- simulate_process(type = "nonlinear")
#' plot(dat$x, dat$y)
#'
#' @export
simulate_process <- function(x = NULL, n_basis = 50, sd_fun = 1, sd = 0.1, sd_poly = 0.1, type = c("linear", "nonlinear", "quadratic", "nondynamic"), p = 1, pred_step = 16, normalize = FALSE){
  type <- match.arg(type)
  if(type == "linear"){
    g <- simulate_linear_function(sd_poly = sd_poly)
  }
  else if(type == "quadratic"){
    g <- simulate_quadratic_function(sd_poly = sd_poly)
  }
  else if(type == "nonlinear") {
    g <- simulate_nonlinear_function(n_basis = n_basis, sd_function = sd_fun, sd_poly = sd_poly, p = p, pred_step = pred_step)
  }
  else if(type == "nondynamic") {
    g <- simulate_nondynamic_function(sd_poly = sd_poly)
  }
  else {
    stop("type must be one of 'linear', 'nonlinear', 'nondynamic'")
  }

  if(normalize){
    if(is.null(x)) {
      x_vec <- seq(0, 16, length.out = 100)
    } else {
      x_vec <- x
    }
    if(diff(range(g(x_vec))) > 2){
      gScale <- function(x) ((g(x) - min(g(x_vec)))/diff(range(g(x_vec))))
      gFinal <- function(x) gScale(x) - mean(gScale(x_vec))
    }
    else{
      gFinal <- function(x) g(x) - mean(g(x_vec))
    }
  }
  else{
    gFinal <- g
  }

  return(simulate_data(g = gFinal, sd = sd, x = x))
}



