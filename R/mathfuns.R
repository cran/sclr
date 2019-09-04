# Shortcut math functions

# Inverse logit function
invlogit <- function(x) exp(x) / (1 + exp(x))

# Commonly occuring exp(Xb) expression
get_exp_Xb <- function(y, x, pars) exp(get_Xb(y, x, pars))

# Result of Xb matrix multiplication
get_Xb <- function(y, x, pars) {
  pars_betas <- get_betas_only(pars)
  xb <- x %*% pars_betas
  return(xb)
}
