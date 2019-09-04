# Testing the predict method and protection calculations 
# (these are powered by predict).
#
# Arseniy Khvorov
# Created 2019/09/02
# Last edit 2019/09/02

library(sclr)

fit1 <- sclr(status ~ logHI, sclr_one_titre_data)
preddata1 <- data.frame(logHI = c(0, 3))
fit2 <- sclr(status ~ logHI + logNI, sclr_two_titre_data)
preddata2 <- data.frame(logHI = c(0, 3), logNI = c(0, 3))
fittrans <- sclr(status ~ log(exp(logHI)), sclr_one_titre_data)

test_that(
  "Predict method works", {
    expect_named(predict(fit1, preddata1))
    expect_named(predict(fit2, preddata2))
  }
)

test_that(
  "Predict method throws an error when used incorrectly", {
    expect_error(predict(fit2, preddata1)) # Not enough covariates
  }
)

test_that(
  "Predict method works with transformations", {
    expect_named(predict(fittrans, preddata1))
  }
)

test_that(
  "Protection level calculation works", {
    expect_named(get_protection_level(fit1, "logHI"))
    expect_named(
      get_protection_level(fit2, "logHI", data.frame(logNI = c(0, 5)))
    )
  }
)

test_that(
  "Protection level calculation works with transformations", {
    expect_named(get_protection_level(fittrans, "logHI"))
  }
)
