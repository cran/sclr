# Functions to compute the score matrix
# Arseniy Khvorov
# Created 2019/09/11
# Last edit 2019/10/15

#' Score matrix
#'
#' @param y Model response
#' @param x Model matrix
#' @param pars Parameter matrix
#' @param exp_Xb exp(xb) term
#'
#' @noRd
get_scores <- function(y, x, pars, exp_Xb) {
  dl <- get_dl(y, x, pars, exp_Xb)
  beta_part <- get_score_beta_part(y, x, pars, exp_Xb)
  scores <- matrix(c(dl, beta_part), ncol = 1)
  rownames(scores) <- paste0("d_", get_par_names(x))
  scores
}

#' First derivative of log-likelihood with respect to lambda
#'
#' @inheritParams get_scores
#'
#' @noRd
get_dl <- function(y, x, pars, exp_Xb) {
  lambda <- pars[1, ]
  dl <- y / lambda - (1 - y) / (1 + exp_Xb - lambda)
  dl <- sum(dl)
  return(dl)
}

#' First derivative of log-likelihood with respect to betas
#'
#' @inheritParams get_scores
#'
#' @noRd
get_score_beta_part <- function(y, x, pars, exp_Xb) {
  common <- get_score_beta_common(y, x, pars, exp_Xb) # A matrix with 1 col
  beta_score_contr <- x * common[, 1] # Non-matrix multiplication
  colSums(beta_score_contr)
}

#' Common part of the beta derivatives
#'
#' @inheritParams get_scores
#'
#' @noRd
get_score_beta_common <- function(y, x, pars, exp_Xb) {
  (1 - y) * get_dbc_f1(y, x, pars, exp_Xb) - get_dbc_f2(y, x, pars, exp_Xb)
}

#' First fraction term of common part of the beta derivatives
#'
#' @inheritParams get_scores
#'
#' @noRd
get_dbc_f1 <- function(y, x, pars, exp_Xb) {
  lambda <- pars[1, ]
  exp_Xb / (1 + exp_Xb - lambda)
}

#' Second fraction term of common part of the beta derivatives
#'
#' @inheritParams get_scores
#'
#' @noRd
get_dbc_f2 <- function(y, x, pars, exp_Xb) {
  exp_Xb / (1 + exp_Xb)
}
