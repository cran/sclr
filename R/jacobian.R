# Functions to caclulate the log-likelihood second derivative matrix.
# Arseniy Khvorov
# Created 2019/07/31
# Last edit 2019/10/15

#' Second derivative matrix
#'
#' @param y Model response
#' @param x Model matrix
#' @param pars_mat Parameter matrix
#' @param exp_Xb exp(xb) term
#' @param x_coeffs Matrix of pairwise products of elements of x
#'
#' @noRd
get_jacobian <- function(y, x, pars_mat, exp_Xb, x_coeffs) {
  
  lambda_entries <- get_lambda_array(y, x, pars_mat, exp_Xb)
  beta_part <- get_jac_beta_part(y, x, pars_mat, exp_Xb, x_coeffs)
  
  n_x <- ncol(x)
  
  # Convert the beta part to a matrix
  beta_part_mat <- build_symm_mat(beta_part)

  # Combine into full matrix
  n_par <- n_x + 1
  nms <- paste0("d_", get_par_names(x))
  sec_dev_mat <- matrix(
    0, nrow = n_par, ncol = n_par, dimnames = list(nms, nms)
  )
  sec_dev_mat[1, ] <- lambda_entries
  sec_dev_mat[, 1] <- lambda_entries
  sec_dev_mat[2:n_par, 2:n_par] <- beta_part_mat

  sec_dev_mat
}

#' lambda row and column entries
#'
#' @param y Model response
#' @param x Model matrix
#' @param pars_mat Parameter matrix
#' @param exp_Xb exp(xb) term
#'
#' @noRd
get_lambda_array <- function(y, x, pars_mat, exp_Xb) {
  dl2 <- get_dl2(y, x, pars_mat, exp_Xb)
  dldbs <- get_dldbs(y, x, pars_mat, exp_Xb)
  c(dl2, dldbs)
}

#' Second derivative of log-likelihood with respect to lambda
#'
#' @inheritParams get_lambda_array
#'
#' @noRd
get_dl2 <- function(y, x, pars_mat, exp_Xb) {
  lambda <- pars_mat[1, ]
  dl2 <- -y / lambda^2 - (1 - y) / (1 + exp_Xb - lambda)^2
  sum(dl2)
}

#' Second derivative of log-likelihood with respect to lambda and betas
#'
#' @inheritParams get_lambda_array
#'
#' @noRd
get_dldbs <- function(y, x, pars_mat, exp_Xb) {
  common <- get_dldbs_common(y, x, pars_mat, exp_Xb)
  dldbs_contr <- x * common[, 1]
  colSums(dldbs_contr)
}

#' Common part of second derivative of log-likelihood with respect to lambda
#' and betas
#'
#' @inheritParams get_lambda_array
#'
#' @noRd
get_dldbs_common <- function(y, x, pars_mat, exp_Xb) {
  lambda <- pars_mat[1, ]
  (1 - y) * exp_Xb / (1 + exp_Xb - lambda)^2
}

#' Second derivative of log-likelihood with respect to betas
#'
#' @inheritParams get_jacobian
#'
#' @noRd
get_jac_beta_part <- function(y, x, pars_mat, exp_Xb, x_coeffs) {
  
  common <- get_jac_beta_common(y, x, pars_mat, exp_Xb)
  
  # Get the unique triangle of the beta part
  beta_contr <- x_coeffs * common[, 1]
  colSums(beta_contr)
}

#' Common part of second derivative of log-likelihood with respect to betas
#'
#' @inheritParams get_lambda_array
#'
#' @noRd
get_jac_beta_common <- function(y, x, pars_mat, exp_Xb) {
  cf1 <- get_jac_beta_common_f1(y, x, pars_mat, exp_Xb)
  cf2 <- get_jac_beta_common_f2(y, x, pars_mat, exp_Xb)
  cf1 - cf2
}

#' First fraction of common part of second derivative of 
#' log-likelihood with respect to betas
#'
#' @inheritParams get_lambda_array
#'
#' @noRd
get_jac_beta_common_f1 <- function(y, x, pars_mat, exp_Xb) {
  lambda <- pars_mat[1, ]
  (1 - y) * (1 - lambda) * exp_Xb / (1 + exp_Xb - lambda)^2
}

#' Second fraction of common part of second derivative of 
#' log-likelihood with respect to betas
#'
#' @inheritParams get_lambda_array
#'
#' @noRd
get_jac_beta_common_f2 <- function(y, x, pars_mat, exp_Xb) {
  exp_Xb / (1 + exp_Xb)^2
}
