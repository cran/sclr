# Likelihood-related functions
# Arseniy Khvorov
# Created 2019/07/31
# Last edit 2019/09/02

#' Log-likelihood
#' 
#' Computes the log-likelihood of the scaled logit model
#' at a given set of parameter estimates 
#' (or the MLE if \code{pars} is not supplied).
#'
#' @param fit An object returned by \code{\link{sclr}}.
#' @param pars A named vector of parameter values. If \code{NULL}
#' then the estimates from \code{fit} will be used.
#' 
#' @export
sclr_log_likelihood <- function(fit, pars = NULL) {
  
  if (is.null(pars)) pars <- fit$parameters
  
  if (!is.matrix(pars)) {
    pars_vec <- pars
    pars <- as.matrix(pars_vec, ncol = 1)
    rownames(pars) <- names(pars_vec)
  }
  
  lambda <- pars["lambda", ]
  
  y <- fit$y
  x <- fit$x
  
  exp_Xb <- get_exp_Xb(y, x, pars)
  
  l <- y * log(lambda) + (1 - y) * log(1 + exp_Xb - lambda) - log(1 + exp_Xb)
  l <- sum(l)
  
  return(l)
}
