# Fitting functions
# Arseniy Khvorov
# Created 2019/07/31
# Last edit 2019/10/16

#' Fitter function for the scaled logit model
#'
#' Computing engine behind \code{\link{sclr}}.
#'
#' The likelihood maximisation uses the Newton-Raphson algorithm. Initial values
#' are always 1 for the covariate coefficients (and the associated intercept)
#' and the proportion of infected for the baseline risk.
#'
#' The algorithm will pick a new guess and restart under a set of conditions.
#'
#' 1) Algorithm's iteration produces estimate guesses that cannot be used -
#' baseline risk outside of (0, 1) since likelihood is undefined.
#'
#' 2) The second derivative matrix produced by the current estimates is "bad" -
#' positive semi-definite or missing values due to failing large 
#' number calculations.
#'
#' @param y A vector of observations.
#' @param x A design matrix.
#' @param tol Tolerance. Used when \code{n_iter} is \code{NULL}.
#' @param n_iter Number of Newton-Raphson iterations. \code{tol} is ignored when
#'   this is not \code{NULL}.
#' @param max_tol_it Maximum tolerated iterations. If it fails to converge
#'   within this number of iterations, will return with an error.
#' @param n_conv Number of times the algorithm has to converge (to work around
#'   local maxima).
#' @param conventional_names If \code{TRUE}, estimated parameter names will be
#'   (Baseline), (Intercept) and the column names in the model matrix. Otherwise
#'   - lambda, beta_0 and beta_ prefix in front of column names in the model
#'   matrix.
#'   
#' @importFrom rlang abort
#'
#' @export
sclr_fit <- function(y, x, 
                     tol = 10^(-7), 
                     n_iter = NULL, max_tol_it = 10^4, n_conv = 3,
                     conventional_names = FALSE) {
  
  # Parameter matrix with initial values
  pars_mat <- get_init_pars_mat(x, y, conventional_names)

  x_coeffs <- get_x_coeffs(x) # To avoid recalculations

  # Work out the MLEs
  n_iter_cur <- 1
  conv_count <- 0
  lls <- c()
  rets <- list()
  while (TRUE) {
    
    # Check that the current set is workable
    if (is_bad_pars(pars_mat)) pars_mat <- guess_again(pars_mat)
    pars_mat_prev <- pars_mat # Save for later
    
    # Calculate the commonly occuring expession
    exp_Xb <- get_exp_Xb(y, x, pars_mat)
    
    # Log-likelihood second derivative (negative of information) matrix
    jacobian_mat <- get_jacobian(y, x, pars_mat, exp_Xb, x_coeffs)
    if (is_bad_jac(jacobian_mat)) {
      pars_mat <- guess_again(pars_mat)
      next
    }
    
    # Invert to get the negative of the covariance matrix
    inv_jacobian_mat <- try(solve(jacobian_mat), silent = TRUE)
    if (inherits(inv_jacobian_mat, "try-error") || 
        any(is.na(inv_jacobian_mat))) {
      pars_mat <- guess_again(pars_mat)
      next
    }
    
    # Generalised Newton-Raphson
    scores_mat <- get_scores(y, x, pars_mat, exp_Xb)
    pars_mat <- pars_mat_prev - inv_jacobian_mat %*% scores_mat
    
    # Check convergence
    if (has_converged(pars_mat, pars_mat_prev, tol)) {
      conv_count <- conv_count + 1
      lls[[conv_count]] <- sclr_log_likelihood(x = x, y = y, pars = pars_mat)
      rets[[conv_count]] <- list(
        "invjac" = inv_jacobian_mat, "pars" = pars_mat
      )
      if (conv_count == n_conv) break
      pars_mat <- guess_again(pars_mat)
      next
    } else {
      if (!is.null(n_iter) && (n_iter_cur >= n_iter)) {
        lls[[conv_count + 1]] <- sclr_log_likelihood(
          x = x, y = y, pars = pars_mat
        )
        rets[[conv_count + 1]] <- list(
          "invjac" = inv_jacobian_mat, "pars" = pars_mat
        )
        break
      }
      if (n_iter_cur >= max_tol_it) {
        abort(paste0("did not converge in ", max_tol_it, " iterations"))
      }
      n_iter_cur <- n_iter_cur + 1
      next
    }
  }
  
  if (all(is.nan(lls))) 
    abort(paste0("no reportable results in ", n_iter_cur, " iterations"))
  
  i_best <- which.max(lls)
  pars_mat <- rets[[i_best]][["pars"]]
  inv_jacobian_mat <- rets[[i_best]][["invjac"]]
  
  # Build the return list
  parameters <- as.vector(pars_mat)
  names(parameters) <- rownames(pars_mat)
  
  covariance_mat <- -inv_jacobian_mat
  dimnames(covariance_mat) <- list(names(parameters), names(parameters))
  list(
    parameters = parameters,
    covariance_mat = covariance_mat,
    n_converge = n_iter_cur
  )
}

#' Initial parameter matrix
#'
#' @param x Model matrix.
#' @param y Model response.
#' @param conventional_names Controls parameter names.
#'
#' @noRd
get_init_pars_mat <- function(x, y, conventional_names) {
  n_par <- ncol(x) + 1
  pars_mat <- matrix(rep(1, n_par))
  rownames(pars_mat) <- get_par_names(x, conventional_names)
  pars_mat[1, ] <- mean(y)
  pars_mat
}

#' Check the parameter matrix.
#'
#' Checks if the current parameter guesses are OK for derivative calculations.
#' Returns \code{TRUE} if they are and \code{FALSE} otherwise.
#'
#' @param pars_mat The current matrix of parameter guesses.
#'
#' @noRd
is_bad_pars <- function(pars_mat) {
  lambda <- pars_mat[1, ]
  (lambda < 0) || (lambda > 1) # Likelihood undefined
}

#' Check the second derivative matrix
#'
#' Checks if the current parameter guesses are OK for derivative calculations.
#' Returns \code{TRUE} if they are and \code{FALSE} otherwise.
#'
#' @param jac The current second derivative matrix.
#'
#' @noRd
is_bad_jac <- function(jac) {
  if (any(is.na(jac))) return(TRUE)
  is_convex(jac)
}

#' Create a new guess
#'
#' Creates a matrix with new parameter guesses.
#' The matrix has the same structure as the matrix with previous guesses.
#'
#' @param pars_mat The current matrix of parameter guesses.
#' 
#' @importFrom stats runif rnorm
#'
#' @noRd
guess_again <- function(pars_mat) {
  delta <- matrix(
    c(
      runif(1, min = 0, max = 1),
      rnorm(nrow(pars_mat) - 1, mean = 0, sd = 2)
    ),
    ncol = 1
  )
  rownames(delta) <- rownames(pars_mat)
  return(delta)
}

#' Check for convexity
#'
#' @param jac Second derivative matrix
#' 
#' @importFrom dplyr near
#'
#' @noRd
is_convex <- function(jac) {
  eigenvals <- eigen(jac, only.values = TRUE)$values
  eigenvals[near(eigenvals, 0)] <- 0
  any(eigenvals > 0)
}

#' Check convergence
#'
#' @param pars_mat Current guess
#' @param pars_mat_prev Last guess
#' @param tol Tolerance
#'
#' @noRd
has_converged <- function(pars_mat, pars_mat_prev, tol) {
  deltas <- abs(pars_mat - pars_mat_prev)
  all(deltas < tol)
}
