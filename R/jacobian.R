# Functions to caclulate the log-likelihood second derivative matrix.
# Arseniy Khvorov
# Created 2019/07/31
# Last edit 2019/09/02

# Computes the log-likelihood second derivative matrix
get_jacobian <- function(y, x, pars_mat, exp_Xb) {
  
  lambda_entries <- get_lambda_array(y, x, pars_mat, exp_Xb)
  beta_part <- get_jac_beta_part(y, x, pars_mat, exp_Xb)
  
  n_x <- ncol(x)
  
  # Convert the beta part to a matrix
  beta_part_mat <- build_symm_mat(beta_part)

  # Combine into full matrix
  n_par <- n_x + 1
  nms <- c("d_lambda", paste0("d_", get_par_names(x)))
  sec_dev_mat <- matrix(
    0, nrow = n_par, ncol = n_par, dimnames = list(nms, nms)
  )
  sec_dev_mat[1, ] <- lambda_entries
  sec_dev_mat[, 1] <- lambda_entries
  sec_dev_mat[2:n_par, 2:n_par] <- beta_part_mat

  return(sec_dev_mat)
}

# Computes the lambda row and column entries
get_lambda_array <- function(y, x, pars_mat, exp_Xb) {
  dl2 <- get_dl2(y, x, pars_mat, exp_Xb)
  dldbs <- get_dldbs(y, x, pars_mat, exp_Xb)
  lambda_entries <- c(dl2, dldbs)
  return(lambda_entries)
}

# Computes the second derivative of log-likelihood with respect to lambda
get_dl2 <- function(y, x, pars_mat, exp_Xb) {
  lambda <- pars_mat["lambda", ]
  dl2 <- -y / lambda^2 - (1 - y) / (1 + exp_Xb - lambda)^2
  dl2 <- sum(dl2)
  return(dl2)
}

# Computes the second derivatives of log-likelihood with respect to lambda and
# betas
get_dldbs <- function(y, x, pars_mat, exp_Xb) {
  common <- get_dldbs_common(y, x, pars_mat, exp_Xb)
  dldbs_contr <- x * common[, 1]
  dldbs <- colSums(dldbs_contr)
  return(dldbs)
}

# Computes the common part of the
# second derivatives of log-likelihood with respect to lambda and betas
get_dldbs_common <- function(y, x, pars_mat, exp_Xb) {
  lambda <- pars_mat["lambda", ]
  common <- (1 - y) * exp_Xb / (1 + exp_Xb - lambda)^2
  return(common)
}

# Computes the second derivatives of log likelihood with respect to the betas
get_jac_beta_part <- function(y, x, pars_mat, exp_Xb) {
  
  common <- get_jac_beta_common(y, x, pars_mat, exp_Xb)
  
  # Create the x coefficients
  x_coeffs <- get_x_coeffs(x)
  
  # Get the unique triangle of the beta part
  beta_contr <- x_coeffs * common[, 1]
  beta_part <- colSums(beta_contr)
  
  return(beta_part)
}

# Computes the common part of the second derivatives
# of log likelihood with respect to the betas
get_jac_beta_common <- function(y, x, pars_mat, exp_Xb) {
  cf1 <- get_jac_beta_common_f1(y, x, pars_mat, exp_Xb)
  cf2 <- get_jac_beta_common_f2(y, x, pars_mat, exp_Xb)
  common <- cf1 - cf2
  return(common)
}

# Computes the first fraction of the common part of the second derivatives
# of log likelihood with respect to the betas
get_jac_beta_common_f1 <- function(y, x, pars_mat, exp_Xb) {
  lambda <- pars_mat["lambda", ]
  common_f1 <- (1 - y) * (1 - lambda) * exp_Xb / (1 + exp_Xb - lambda)^2
  return(common_f1)
}

# Computes the second fraction of the common part of the second derivatives
# of log likelihood with respect to the betas
get_jac_beta_common_f2 <- function(y, x, pars_mat, exp_Xb) {
  common_f2 <- exp_Xb / (1 + exp_Xb)^2
  return(common_f2)
}
