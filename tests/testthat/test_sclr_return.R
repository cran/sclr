# Tests of sclr returning the expected values.
# 
# These 'correct' values are coming from the verion of sclr that was
# verified to be good through simulations.
#
# Arseniy Khvorov
# Created 2019/09/02
# Last edit 2019/09/02

library(sclr)

fit1 <- sclr(status ~ logHI, sclr_one_titre_data)
expected_ll1 <- -2292.91615285242
expected_est1 <- c(
  "lambda" = 0.24382797620260,
  "beta_0" = -7.76395203110848,
  "beta_logHI" = 2.08804817370896
)

fit2 <- sclr(status ~ logHI + logNI, sclr_two_titre_data)
expected_ll2 <- -774.953806738847
expected_est2 <- c(
  "lambda" = 0.240734189154097,
  "beta_0" = -8.483447969938824,
  "beta_logHI" = 2.347614894467592,
  "beta_logNI" = 2.181710534930221
)

test_that(
  "Returns the expected maximum likelihood", {
    ll1 <- sclr_log_likelihood(fit1)
    expect_equal(ll1, expected_ll1)
    ll2 <- sclr_log_likelihood(fit2)
    expect_equal(ll2, expected_ll2)
  }
)

test_that(
  "Returns the expected parameter estimates", {
    est1 <- coef(fit1)
    expect_equal(est1, expected_est1)
    est2 <- coef(fit2)
    expect_equal(est2, expected_est2)
  }
)
