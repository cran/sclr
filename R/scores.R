# Functions to compute the score matrix

# Computes the score matrix
get_scores <- function(y, x, pars, exp_Xb) {
  dl <- get_dl(y, x, pars, exp_Xb)
  beta_part <- get_score_beta_part(y, x, pars, exp_Xb)
  scores <- matrix(c(dl, beta_part), ncol = 1)
  rownames(scores) <- c("d_lambda", paste0("d_", get_par_names(x)))
  return(scores)
}

# Computes the first derivative of log-likelihood with respect to lambda
get_dl <- function(y, x, pars, exp_Xb) {
  lambda <- pars["lambda", ]
  dl <- y / lambda - (1 - y) / (1 + exp_Xb - lambda)
  dl <- sum(dl)
  return(dl)
}

# Computes the first derivatives of log-likelihood with respect to the betas
get_score_beta_part <- function(y, x, pars, exp_Xb) {
  common <- get_score_beta_common(y, x, pars, exp_Xb) # A matrix with 1 col
  beta_score_contr <- x * common[, 1] # Non-matrix multiplication
  beta_score <- colSums(beta_score_contr)
  return(beta_score)
}

# Computes common part of the beta derivatives
get_score_beta_common <- function(y, x, pars, exp_Xb) {
  bc <- (1 - y) * get_dbc_f1(y, x, pars, exp_Xb) - 
    get_dbc_f2(y, x, pars, exp_Xb)
  return(bc)
}

# Computes the first fraction term of common
get_dbc_f1 <- function(y, x, pars, exp_Xb) {
  lambda <- pars["lambda", ]
  f1 <- exp_Xb / (1 + exp_Xb - lambda)
  return(f1)
}

# Computes the second fraction term of common
get_dbc_f2 <- function(y, x, pars, exp_Xb) {
  f2 <- exp_Xb / (1 + exp_Xb)
  return(f2)
}
