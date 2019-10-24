# sclr class functions
# Arseniy Khvorov
# Created 2019/07/31
# Last edit 2019/10/15

#' Fits the scaled logit model
#'
#' Used to fit the scaled logit model from Dunning (2006).
#'
#' The model is of the form \deqn{P(Y = 1) = \lambda(1 - logit^{-1}(\beta_0 +
#' \beta_1X_1 + \beta_2X_2 + ... + \beta_kX_k))} Where \eqn{Y} is the binary
#' outcome indicator, (eg. 1 - infected, 0 - not infected). \eqn{X} - covariate.
#' \eqn{k} - number of covariates. Computing engine behind the fitting is
#' \code{\link{sclr_fit}}.
#'
#' @param formula an object of class "formula": a symbolic description of the
#'   model to be fitted.
#' @param data a data frame.
#' @param ci_lvl Confidence interval level for the parameter estimates.
#' @inheritParams sclr_fit
#'
#' @return An object of class \code{sclr}. This is a list with the following
#'   elements:
#'
#'   \item{parameters}{Maximum likelihood estimates of the parameter values.}
#'
#'   \item{covariance_mat}{The variance-covariance matrix of the parameter
#'   estimates.}
#'
#'   \item{n_converge}{The number of Newton-Raphson iterations (including
#'   resets) that were required for convergence.}
#'
#'   \item{x}{Model matrix derived from \code{formula} and \code{data}.}
#'
#'   \item{y}{Response matrix derived from \code{formula} and \code{data}.}
#'   
#'   \item{call}{The original call to \code{sclr}.}
#'
#'   \item{model}{Model frame object derived from \code{formula} and 
#'   \code{data}.}
#'
#'   \item{terms}{Terms object derived from model frame.}
#'   
#'   \item{ci}{Confidence intervals of the parameter estimates.}
#'   
#'   \item{log_likelihood}{Value of log-likelihood calculated at the ML
#'   estimates of parameters.} 
#'
#'   Methods supported: \code{\link[=print.sclr]{print}},
#'   \code{\link[=vcov.sclr]{vcov}}, \code{\link[=coef.sclr]{coef}},
#'   \code{\link[=model.frame.sclr]{model.frame}},
#'   \code{\link[=model.matrix.sclr]{model.matrix}},
#'   \code{\link[=summary.sclr]{summary}}, \code{\link[=predict.sclr]{predict}},
#'   \code{\link[=tidy.sclr]{tidy}} (\code{\link{broom}} package).
#'
#' @references Dunning AJ (2006). "A model for immunological correlates of
#'   protection." Statistics in Medicine, 25(9), 1485-1497.
#'   \url{https://doi.org/10.1002/sim.2282}.
#'
#' @examples
#' library(sclr)
#' fit1 <- sclr(status ~ logHI, one_titre_data)
#' summary(fit1)
#' @importFrom stats confint model.frame model.matrix model.response
#' @importFrom rlang abort
#'
#' @export
sclr <- function(formula, data = NULL, 
                 ci_lvl = 0.95, 
                 tol = 10^(-7), 
                 n_iter = NULL, 
                 max_tol_it = 10^4, 
                 n_conv = 3,
                 conventional_names = FALSE) {
  
  if (!inherits(formula, "formula") || missing(formula)) 
    abort("must supply a formula")

  cl <- match.call()
  mf <- model.frame(formula, data, drop.unused.levels = TRUE)

  # Design matirix
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf)

  # Response vector
  y <- model.response(mf)
  if (!all(y %in% c(0, 1))) abort("response should be a vector with 0 and 1")
  if (!is.numeric(y)) abort("response should be numeric")

  # Actual model fit
  fit <- sclr_fit(y, x, tol, n_iter, max_tol_it, n_conv, conventional_names)

  # Build the return list
  fit <- new_sclr(fit, x, y, cl, mf, mt)
  fit$ci <- confint(fit, level = ci_lvl)
  fit$log_likelihood <- sclr_log_likelihood(fit)
  fit
}

#' Create a new \code{sclr} object
#' 
#' \code{new_sclr} creates the object \code{\link{sclr}} returns.
#' \code{is_sclr} checks if the object is of class \code{sclr}.
#'
#' @param fit A list returned by \code{\link{sclr_fit}}.
#' @param x Model matrix.
#' @param y Model response.
#' @param cl Call.
#' @param mf Model frame.
#' @param mt Model terms.
#'
#' @return \code{sclr} object
#' @export
new_sclr <- function(fit, x, y, cl, mf, mt) {
  stopifnot(is.list(fit))
  stopifnot(is.matrix(x))
  stopifnot(is.numeric(y))
  stopifnot(is.call(cl))
  stopifnot(is.data.frame(mf))
  stopifnot(inherits(mt, "terms"))
  fit <- c(fit, list(
    x = x, y = y,
    call = cl, model = mf, terms = mt
  ))
  class(fit) <- "sclr"
  fit
}

#' @rdname new_sclr
#' @export
is_sclr <- function(fit) "sclr" %in% class(fit)
