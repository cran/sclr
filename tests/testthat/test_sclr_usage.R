# Tests of sclr returning/throwing errors properly
# Arseniy Khvorov
# Created 2019/09/02
# Last edit 2019/09/02

library(sclr)

test_that(
  "Returns a named object when supposed to", {
    expect_named(sclr(status ~ logHI, sclr_one_titre_data))
    expect_named(sclr(status ~ logHI + logNI, sclr_two_titre_data))
  }
)

test_that(
  "Error with missing parameters", {
    expect_error(sclr(status ~ logHI))
    expect_error(sclr(sclr_one_titre_data))
  }
)

test_that(
  "Error with unexpected outcome", {
    # Factor
    dat <- sclr_one_titre_data
    dat$status <- as.factor(dat$status)
    expect_error(sclr(status ~ logHI, dat))
    # 1's and 2's
    dat$status <- as.numeric(dat$status)
    expect_error(sclr(status ~ logHI, dat))
  }
)
