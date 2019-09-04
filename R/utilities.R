# General utility functions

# Returns parameter names based on the model matrix
get_par_names <- function(x) {
  par_names <- c("beta_0", paste0("beta_", colnames(x)[-1]))
  return(par_names)
}

# Narrows the parameter list to betas
get_betas_only <- function(pars_mat) {
  pars_betas <- pars_mat[grepl("beta", rownames(pars_mat)), ]
  pars_betas_mat <- matrix(pars_betas, ncol = 1)
  rownames(pars_betas_mat) <- names(pars_betas)
  return(pars_betas_mat)
}

# Creates a matrix of all products of pairwise multiplication of elements of x
get_x_coeffs <- function(x) {
  n_x <- ncol(x)
  n_coefs <- choose(n_x, 2) + n_x
  x_coeffs <- matrix(0, ncol = n_coefs, nrow = nrow(x))
  cur_coeff <- 1
  for (i in 1:n_x) {
    for (j in 1:n_x) {
      if (j < i) next
      x_coeffs[, cur_coeff] <- x[, i] * x[, j]
      cur_coeff <- cur_coeff + 1
    }
  }
  return(x_coeffs)
}

# Creates a symmeterical matrix from vectors of unique triangles
build_symm_mat <- function(obj) {
  if (is.vector(obj)) return(build_symm_mat_one(obj))
  if (!is.data.frame(obj) & !is.matrix(obj)) 
    stop("obj needs to be a vector, a matrix or a dataframe")
  obj <- split(obj, seq(nrow(obj)))
  mats <- lapply(obj, build_symm_mat_one)
  return(mats)
}

# Creates a symmeterical matrix from a vector of the unique triangle
build_symm_mat_one <- function(triang_vec) {
  matdim <- get_symm_dims(length(triang_vec))
  mat <- matrix(0, nrow = matdim, ncol = matdim)
  mat[lower.tri(mat, diag = TRUE)] <- triang_vec
  mat <- t(mat)
  mat[lower.tri(mat, diag = TRUE)] <- triang_vec
  return(mat)
}

# Returns the dimensions of the symmetrical matrix given the number of
# elements in its triangle
get_symm_dims <- function(len) {
  rtexp <- 1 + 8 * len
  rt <- as.integer(sqrt(rtexp))
  if ((rt^2 - rtexp) > 10^(-7)) stop(len, " is not a valid length of triangle")
  return((-1 + rt) / 2)
}

# Returns the message to indicate mismatch of expected and provided variables
nvar_mismatch_msg <- function(name, nvar_exp, nvar_got) {
  paste("in", name, "expected", nvar_exp, "variables", "got", nvar_got)
}
