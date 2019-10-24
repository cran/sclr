# sclr 0.1.0-dev

- Added `sclr_ideal_data` function to simulate ideal data for the model.

- Made simulations in data-raw self-contained.

- Added the ability to return parameter names that are more conventional (e.g.
"(Intercept)" instead of "beta_0"). See `conventional_names` argument in
`?sclr`.

- Made convergence stricter to avoid local maxima. Argument `n_conv` to `sclr`
and `sclr_fit` sets the number of times the algorithm has to converge. Best set
(the one with maximum likelihood) is chosen out of `n_conv` sets. Previously,
the algorithm only converged once.

- `sclr_log_likelihood` can now be called with a model matrix and a model response.

- Minor performance optimisations.

# sclr 0.1.0

First release.

## Main features

- Fits the scaled logit model using the Newton-Raphson method.

- Supports the predict method for the expected value of the linear beta X part
of the model.

- Can look for covariate values corresponding to a particular protection level.