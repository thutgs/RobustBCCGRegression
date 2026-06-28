test_that("mdpde_bccg with alpha=0 returns an mdpde_bccg object", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  set.seed(123)
  n   <- 60
  x   <- runif(n)
  y   <- gamlss.dist::rBCCG(n, mu = exp(3 + x), sigma = 0.1, nu = 1.5)
  dat <- data.frame(y = y, x = x)

  fit <- suppressMessages(mdpde_bccg(y ~ x, data = dat, alpha = 0))
  expect_s3_class(fit, "mdpde_bccg")
  expect_equal(fit$alpha, 0)
  expect_length(fit$coefficients$mu, 2)   # intercept + x
  expect_length(fit$coefficients$sigma, 1)
  expect_length(fit$coefficients$nu, 1)
})

test_that("mdpde_bccg with alpha>0 returns a converged mdpde_bccg object", {
  skip_if_not_installed("gamlss.dist")
  skip_if_not_installed("robustbase")

  set.seed(42)
  n   <- 80
  x   <- runif(n)
  y   <- gamlss.dist::rBCCG(n, mu = exp(3 + x), sigma = 0.1, nu = 2)
  dat <- data.frame(y = y, x = x)

  fit <- suppressMessages(mdpde_bccg(y ~ x, data = dat, alpha = 0.1))
  expect_s3_class(fit, "mdpde_bccg")
  expect_equal(fit$convergence, 0)
  expect_true(all(is.finite(coef(fit))))
})

test_that("coef.mdpde_bccg respects the parameter argument", {
  skip_if_not_installed("gamlss.dist")

  set.seed(7)
  n   <- 60
  x   <- runif(n)
  y   <- gamlss.dist::rBCCG(n, mu = exp(3 + x), sigma = 0.1, nu = 1)
  dat <- data.frame(y = y, x = x)

  fit <- suppressMessages(mdpde_bccg(y ~ x, data = dat, alpha = 0))
  expect_length(coef(fit, "mu"),    2)
  expect_length(coef(fit, "sigma"), 1)
  expect_length(coef(fit, "nu"),    1)
  expect_length(coef(fit, "all"),   4)
})

test_that("predict.mdpde_bccg returns correct output for each type", {
  skip_if_not_installed("gamlss.dist")

  set.seed(10)
  n   <- 60
  x   <- runif(n)
  y   <- gamlss.dist::rBCCG(n, mu = exp(3 + x), sigma = 0.1, nu = 1)
  dat <- data.frame(y = y, x = x)

  fit <- suppressMessages(mdpde_bccg(y ~ x, data = dat, alpha = 0))

  pred_resp <- predict(fit, type = "response")
  expect_s3_class(pred_resp, "data.frame")
  expect_true(all(c("fit", "se.fit", "lwr", "upr") %in% names(pred_resp)))

  pred_par <- predict(fit, type = "parameter")
  expect_true(all(c("mu", "sigma", "nu") %in% names(pred_par)))

  pred_q <- predict(fit, type = "quantile")
  expect_true(all(c("fit", "lwr_pred", "upr_pred") %in% names(pred_q)))
})

test_that("quantile residuals have approximately zero mean and unit variance", {
  skip_if_not_installed("gamlss.dist")

  set.seed(99)
  n   <- 200
  x   <- runif(n)
  y   <- gamlss.dist::rBCCG(n, mu = exp(3 + x), sigma = 0.1, nu = 1.5)
  dat <- data.frame(y = y, x = x)

  fit <- suppressMessages(mdpde_bccg(y ~ x, data = dat, alpha = 0))
  r   <- fit$residuals
  expect_lt(abs(mean(r)),   0.3)
  expect_lt(abs(var(r) - 1), 0.5)
})
