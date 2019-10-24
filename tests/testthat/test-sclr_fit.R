# Tests of fit-related functions
# Arseniy Khvorov
# Created 2019/10/16
# Last edit 2019/10/16

library(sclr)

mf <- model.frame(status ~ logHI, one_titre_data)
x <- model.matrix(mf, data = one_titre_data)
y <- model.response(mf)

test_that("sclr_fit can be used directly", {
  expect_named(sclr_fit(y, x), c("parameters", "covariance_mat", "n_converge"))
})

test_that("Parameter matrix initalisation and resetting works", {
  init_mat <- get_init_pars_mat(x, y, conventional_names = FALSE)
  expected <- as.matrix(c(mean(y), 1, 1))
  rownames(expected) <- c("lambda", "beta_0", "beta_logHI")
  expect_equal(init_mat, expected)
  expect_true(all(guess_again(init_mat) != guess_again(init_mat)))
})

test_that("Bad guesses are detected", {
  expect_true(is_bad_pars(as.matrix(c(-1, 0, 0))))
  expect_true(is_bad_pars(as.matrix(c(2, 0, 0))))
  expect_true(!is_bad_pars(as.matrix(c(0.5, 0, 0))))
})

small_sample <- sclr_ideal_data(n = 50, seed = 20191024)

test_that("Does not converge with few allowed iterations", {
  expect_error(sclr(status ~ logHI, small_sample, max_tol_it = 10))
})

test_that("Reports last results with few iterations", {
  set.seed(20191024)
  expect_true(is_sclr(sclr(status ~ logHI, small_sample, n_iter = 10)))
})

test_that("Error if no reportable results with few iterations", {
  set.seed(7)
  expect_error(sclr(status ~ logHI, small_sample, n_iter = 10))
})
