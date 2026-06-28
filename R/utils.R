# ==============================================================================
# utils.R
# Internal helper functions: link functions, z/w transforms, score components,
# second derivatives, observed information matrix, and MDPDE gradient.
# ==============================================================================

# ------------------------------------------------------------------------------
# Link functions
# ------------------------------------------------------------------------------

#' Inverse link function
#' @param eta Linear predictor vector.
#' @param link Character string: `"log"`, `"identity"`, or `"sqrt"`.
#' @return Transformed parameter vector.
#' @keywords internal
link_inverse <- function(eta, link) {
  if (link == "log")      return(exp(eta))
  if (link == "identity") return(eta)
  if (link == "sqrt")     return(eta^2)
  stop("Link not supported: ", link)
}

#' First derivative of the inverse link (d mu / d eta)
#' @param theta Parameter vector (on the response scale).
#' @param link Character string.
#' @return Derivative vector.
#' @keywords internal
link_derivative <- function(theta, link) {
  if (link == "log")      return(theta)
  if (link == "identity") return(rep(1, length(theta)))
  if (link == "sqrt")     return(2 * sqrt(theta))
  stop("Link not supported: ", link)
}

#' Second derivative of the inverse link (d^2 mu / d eta^2)
#' @param theta Parameter vector (on the response scale).
#' @param link Character string.
#' @return Second derivative vector.
#' @keywords internal
link_second_derivative <- function(theta, link) {
  if (link == "log")      return(theta)
  if (link == "identity") return(rep(0, length(theta)))
  if (link == "sqrt")     return(rep(2, length(theta)))
  stop("Link not supported: ", link)
}

# ------------------------------------------------------------------------------
# BCCG building blocks: z and w
# ------------------------------------------------------------------------------

#' Standardised residual z for the BCCG distribution
#'
#' Computes z = ((y/mu)^nu - 1) / (nu * sigma) for nu != 0,
#' and z = log(y/mu) / sigma for nu == 0.
#'
#' @param y  Observed response (scalar).
#' @param mu Median parameter (scalar, > 0).
#' @param sigma Approximate CV parameter (scalar, > 0).
#' @param nu Skewness parameter (scalar).
#' @return Scalar standardised residual.
#' @keywords internal
calc_z <- function(y, mu, sigma, nu) {
  if (abs(nu) < 1e-10) {
    return(log(y / mu) / sigma)
  }
  return(((y / mu)^nu - 1) / (nu * sigma))
}

#' Truncation correction ratio phi(w) / Phi(w)
#'
#' Used in the score and Hessian components for sigma and nu when nu != 0.
#'
#' @param sigma Approximate CV parameter.
#' @param nu Skewness parameter (nu != 0).
#' @return Scalar ratio dnorm(1/(sigma|nu|)) / pnorm(1/(sigma|nu|)).
#' @keywords internal
calc_w <- function(sigma, nu) {
  t <- 1 / (sigma * abs(nu))
  return(dnorm(t) / pnorm(t))
}

# ------------------------------------------------------------------------------
# Score functions (individual observation, distribution-level)
# ------------------------------------------------------------------------------

#' Score function for mu (individual observation)
#' @param y,mu,sigma,nu BCCG parameters (scalars).
#' @return Scalar dl/dmu.
#' @keywords internal
score_mu <- function(y, mu, sigma, nu) {
  z <- calc_z(y, mu, sigma, nu)
  if (abs(nu) < 1e-10) return(z / (sigma * mu))
  return((1 / mu) * (nu * z^2 + z / sigma - nu))
}

#' Score function for sigma (individual observation)
#' @param y,mu,sigma,nu BCCG parameters (scalars).
#' @return Scalar dl/dsigma.
#' @keywords internal
score_sigma <- function(y, mu, sigma, nu) {
  z <- calc_z(y, mu, sigma, nu)
  if (abs(nu) < 1e-10) return((z^2 - 1) / sigma)
  w <- calc_w(sigma, nu)
  return((z^2 - 1) / sigma + 1 / (sigma^2 * abs(nu)) * w)
}

#' Score function for nu (individual observation)
#'
#' Returns `NA_real_` when nu == 0 (log-normal case) and is vectorised
#' so that `score_nu(y_vec, mu_vec, sigma_vec, nu_vec)` works.
#'
#' @param y,mu,sigma,nu BCCG parameters (scalars).
#' @return Scalar dl/dnu, or NA when nu == 0.
#' @keywords internal
score_nu <- function(y, mu, sigma, nu) {
  if (abs(nu) < 1e-10) return(NA_real_)
  if (y == 0)          return(-Inf)
  if (is.infinite(y))  return(-sign(nu) * Inf)

  z     <- calc_z(y, mu, sigma, nu)
  w     <- calc_w(sigma, nu)
  term1 <- log(y / mu)
  term2 <- (z / (nu^2 * sigma)) * ((1 + nu * sigma * z) * (nu * log(y / mu) - 1) + 1)
  term3 <- sign(nu) / (sigma * nu^2) * w
  return(term1 - term2 + term3)
}
score_nu <- Vectorize(score_nu)

#' Full score vector for a single observation (regression level)
#'
#' Chains distribution-level scores through the link derivatives and the
#' design-matrix rows to produce the contribution to the gradient of the
#' log-likelihood (or MDPDE objective) for observation i.
#'
#' @param y_i,mu_i,sigma_i,nu_i BCCG parameters for observation i.
#' @param x_mu_i,x_sigma_i,x_nu_i Design-matrix rows (numeric vectors).
#' @param mu_link,sigma_link,nu_link Link function strings.
#' @return Numeric vector of length p_mu + p_sigma + p_nu.
#' @keywords internal
calculate_score <- function(y_i, mu_i, sigma_i, nu_i,
                            x_mu_i, x_sigma_i, x_nu_i,
                            mu_link, sigma_link, nu_link) {
  u_mu    <- score_mu(y_i, mu_i, sigma_i, nu_i)
  u_sigma <- score_sigma(y_i, mu_i, sigma_i, nu_i)
  u_nu    <- if (abs(nu_i) > 1e-6) score_nu(y_i, mu_i, sigma_i, nu_i) else 0

  dmu_deta    <- link_derivative(mu_i,    mu_link)
  dsigma_deta <- link_derivative(sigma_i, sigma_link)
  dnu_deta    <- link_derivative(nu_i,    nu_link)

  grad_mu    <- u_mu    * dmu_deta    * x_mu_i
  grad_sigma <- u_sigma * dsigma_deta * x_sigma_i
  grad_nu    <- u_nu    * dnu_deta    * x_nu_i

  return(c(grad_mu, grad_sigma, grad_nu))
}

# ------------------------------------------------------------------------------
# Second derivatives of the log-likelihood (distribution level)
# ------------------------------------------------------------------------------

#' @keywords internal
ell_mumu <- function(y, mu, sigma, nu) {
  z <- calc_z(y, mu, sigma, nu)
  if (abs(nu) < 1e-10) return((-z / (sigma * mu^2)) - (1 / (sigma^2 * mu^2)))
  R     <- 1 + nu * sigma * z
  term1 <- nu * z^2 + (z / sigma) - nu
  term2 <- (R / sigma) * (2 * nu * z + (1 / sigma))
  return(-(term1 + term2) / mu^2)
}

#' @keywords internal
ell_musigma <- function(y, mu, sigma, nu) {
  z <- calc_z(y, mu, sigma, nu)
  if (abs(nu) < 1e-10) return(-2 * z / (sigma^2 * mu))
  return(-(2 / (mu * sigma)) * (nu * z^2 + (z / sigma)))
}

#' @keywords internal
ell_munu <- function(y, mu, sigma, nu) {
  if (abs(nu) < 1e-10) return(NA_real_)
  z  <- calc_z(y, mu, sigma, nu)
  R  <- 1 + nu * sigma * z
  L  <- log(R)
  dz_dnu <- (R * L - nu * sigma * z) / (nu^2 * sigma)
  return((1 / mu) * (z^2 - 1 + (2 * nu * z + (1 / sigma)) * dz_dnu))
}

#' @keywords internal
ell_sigmasigma <- function(y, mu, sigma, nu) {
  z <- calc_z(y, mu, sigma, nu)
  if (abs(nu) < 1e-10) return(-(3 * z^2 - 1) / sigma^2)
  w  <- calc_w(sigma, nu)
  t  <- 1 / (sigma * abs(nu))
  term1 <- -(3 * z^2 - 1) / sigma^2
  term2 <- -2 * w / (sigma^3 * abs(nu))
  term3 <- w * (t + w) / (sigma^4 * nu^2)
  return(term1 + term2 + term3)
}

#' @keywords internal
ell_sigmanu <- function(y, mu, sigma, nu) {
  if (abs(nu) < 1e-10) return(NA_real_)
  z  <- calc_z(y, mu, sigma, nu)
  R  <- 1 + nu * sigma * z
  L  <- log(R)
  w  <- calc_w(sigma, nu)
  t  <- 1 / (sigma * abs(nu))
  dz_dnu <- (R * L - nu * sigma * z) / (nu^2 * sigma)
  term1 <- (2 * z / sigma) * dz_dnu
  term2 <- -sign(nu) * w / (sigma^2 * nu^2)
  term3 <- w * (t + w) / (sigma^3 * nu^3)
  return(term1 + term2 + term3)
}

#' @keywords internal
ell_nunu <- function(y, mu, sigma, nu) {
  if (abs(nu) < 1e-10) return(NA_real_)
  z  <- calc_z(y, mu, sigma, nu)
  R  <- 1 + nu * sigma * z
  L  <- log(R)
  G  <- R * (L - 1) + 1
  w  <- calc_w(sigma, nu)
  t  <- 1 / (sigma * abs(nu))
  term1 <- -G^2 / (nu^4 * sigma^2)
  term2 <- -z * (R * (L^2 - 2 * L + 2) - 2) / (nu^3 * sigma)
  term3 <- -2 * sign(nu) * w / (sigma * nu^3)
  term4 <- w * (t + w) / (sigma^2 * nu^4)
  return(term1 + term2 + term3 + term4)
}

#' Observed Fisher information matrix for a single observation
#' @keywords internal
observed_information_matrix <- function(y_i, mu_i, sigma_i, nu_i,
                                        x_mu_i, x_sigma_i, x_nu_i,
                                        mu_link, sigma_link, nu_link) {
  ell11 <- ell_mumu(y_i, mu_i, sigma_i, nu_i)
  ell12 <- ell_musigma(y_i, mu_i, sigma_i, nu_i)
  ell13 <- ell_munu(y_i, mu_i, sigma_i, nu_i)
  ell22 <- ell_sigmasigma(y_i, mu_i, sigma_i, nu_i)
  ell23 <- ell_sigmanu(y_i, mu_i, sigma_i, nu_i)
  ell33 <- ell_nunu(y_i, mu_i, sigma_i, nu_i)

  u_mu    <- score_mu(y_i, mu_i, sigma_i, nu_i)
  u_sigma <- score_sigma(y_i, mu_i, sigma_i, nu_i)
  u_nu    <- if (abs(nu_i) > 1e-6) score_nu(y_i, mu_i, sigma_i, nu_i) else 0

  w1 <- link_derivative(mu_i,    mu_link)
  w2 <- link_derivative(sigma_i, sigma_link)
  w3 <- link_derivative(nu_i,    nu_link)

  w1p <- link_second_derivative(mu_i,    mu_link)
  w2p <- link_second_derivative(sigma_i, sigma_link)
  w3p <- link_second_derivative(nu_i,    nu_link)

  a11 <- w1^2 * ell11 - u_mu    * w1p
  a22 <- w2^2 * ell22 - u_sigma * w2p
  a33 <- w3^2 * ell33 - u_nu    * w3p
  a12 <- w1 * w2 * ell12
  a13 <- w1 * w3 * ell13
  a23 <- w2 * w3 * ell23

  p1 <- length(x_mu_i); p2 <- length(x_sigma_i); p3 <- length(x_nu_i)
  p  <- p1 + p2 + p3

  x1 <- matrix(x_mu_i,    ncol = 1)
  x2 <- matrix(x_sigma_i, ncol = 1)
  x3 <- matrix(x_nu_i,    ncol = 1)

  H <- matrix(0, p, p)
  H[1:p1,         1:p1]         <- a11 * (x1 %*% t(x1))
  H[1:p1,         (p1+1):(p1+p2)]       <- a12 * (x1 %*% t(x2))
  H[1:p1,         (p1+p2+1):p]  <- a13 * (x1 %*% t(x3))
  H[(p1+1):(p1+p2), 1:p1]               <- a12 * (x2 %*% t(x1))
  H[(p1+1):(p1+p2), (p1+1):(p1+p2)]    <- a22 * (x2 %*% t(x2))
  H[(p1+1):(p1+p2), (p1+p2+1):p]       <- a23 * (x2 %*% t(x3))
  H[(p1+p2+1):p,  1:p1]                 <- a13 * (x3 %*% t(x1))
  H[(p1+p2+1):p,  (p1+1):(p1+p2)]      <- a23 * (x3 %*% t(x2))
  H[(p1+p2+1):p,  (p1+p2+1):p]         <- a33 * (x3 %*% t(x3))

  return(-H)
}

# ------------------------------------------------------------------------------
# MDPDE gradient (analytic, passed to optim)
# ------------------------------------------------------------------------------

#' Analytic gradient of the MDPDE objective for BCCG regression
#'
#' Used internally by [mdpde_bccg()] as the `gr` argument to [stats::optim()].
#'
#' @param theta Numeric parameter vector (betas concatenated).
#' @param y Response vector.
#' @param X_mu,X_sigma,X_nu Design matrices.
#' @param alpha Tuning parameter (> 0).
#' @param mu_link,sigma_link,nu_link Link strings.
#' @param subdivisions Passed to [stats::integrate()].
#' @return Numeric gradient vector of length p_total.
#' @keywords internal
mdpde_gradient <- function(theta, y, X_mu, X_sigma, X_nu, alpha,
                           mu_link, sigma_link, nu_link,
                           subdivisions = 500) {
  p_mu    <- ncol(X_mu); p_sigma <- ncol(X_sigma); p_nu <- ncol(X_nu)
  n       <- length(y)
  p_total <- p_mu + p_sigma + p_nu

  beta_mu    <- theta[1:p_mu]
  beta_sigma <- theta[(p_mu + 1):(p_mu + p_sigma)]
  beta_nu    <- theta[(p_mu + p_sigma + 1):p_total]

  mu_i    <- link_inverse(as.vector(X_mu    %*% beta_mu),    mu_link)
  sigma_i <- link_inverse(as.vector(X_sigma %*% beta_sigma), sigma_link)
  nu_i    <- link_inverse(as.vector(X_nu    %*% beta_nu),    nu_link)

  if (any(!is.finite(mu_i)) || any(!is.finite(sigma_i)) ||
      any(!is.finite(nu_i)) || any(sigma_i <= 1e-6) || any(mu_i <= 1e-6)) {
    return(rep(1e10, p_total))
  }

  d_mu    <- link_derivative(mu_i,    mu_link)
  d_sigma <- link_derivative(sigma_i, sigma_link)
  d_nu    <- link_derivative(nu_i,    nu_link)

  d_obs <- suppressWarnings(
    gamlss.dist::dBCCG(y, mu = mu_i, sigma = sigma_i, nu = nu_i)
  )
  d_obs[!is.finite(d_obs)] <- 0
  d_obs <- pmax(d_obs, 1e-100)
  w_obs <- d_obs^alpha

  u_mu_obs <- u_sigma_obs <- u_nu_obs <- numeric(n)
  for (i in seq_len(n)) {
    u_mu_obs[i]    <- suppressWarnings(score_mu(y[i],    mu_i[i], sigma_i[i], nu_i[i]))
    u_sigma_obs[i] <- suppressWarnings(score_sigma(y[i], mu_i[i], sigma_i[i], nu_i[i]))
    u_nu_obs[i]    <- suppressWarnings(score_nu(y[i],    mu_i[i], sigma_i[i], nu_i[i]))
    if (!is.finite(u_mu_obs[i]))    u_mu_obs[i]    <- 0
    if (!is.finite(u_sigma_obs[i])) u_sigma_obs[i] <- 0
    if (!is.finite(u_nu_obs[i]))    u_nu_obs[i]    <- 0
  }

  emp_mu    <- u_mu_obs    * w_obs
  emp_sigma <- u_sigma_obs * w_obs
  emp_nu    <- u_nu_obs    * w_obs

  grad_eta_mu <- grad_eta_sigma <- grad_eta_nu <- numeric(n)
  const <- (1 + alpha)

  lower_global <- max(1e-6, min(y) * 0.5)
  upper_global <- max(y) * 2

  for (i in seq_len(n)) {
    joint_integrand <- function(t) {
      d_t <- suppressWarnings(
        gamlss.dist::dBCCG(t, mu = mu_i[i], sigma = sigma_i[i], nu = nu_i[i])
      )
      if (!is.finite(d_t) || d_t < 1e-50) return(c(0, 0, 0))
      d_pow    <- d_t^(1 + alpha)
      u_mu_t    <- suppressWarnings(score_mu(t,    mu_i[i], sigma_i[i], nu_i[i]))
      u_sigma_t <- suppressWarnings(score_sigma(t, mu_i[i], sigma_i[i], nu_i[i]))
      u_nu_t    <- suppressWarnings(score_nu(t,    mu_i[i], sigma_i[i], nu_i[i]))
      if (!is.finite(u_mu_t))    u_mu_t    <- 0
      if (!is.finite(u_sigma_t)) u_sigma_t <- 0
      if (!is.finite(u_nu_t))    u_nu_t    <- 0
      return(c(u_mu_t, u_sigma_t, u_nu_t) * d_pow)
    }

    lower_i <- max(lower_global, mu_i[i] * 0.1)
    upper_i <- min(upper_global, mu_i[i] * 5)

    safe_int <- function(idx) {
      r <- try(
        integrate(function(t) joint_integrand(t)[idx],
                  lower = lower_i, upper = upper_i,
                  subdivisions = subdivisions, rel.tol = 1e-4)$value,
        silent = TRUE
      )
      if (is.numeric(r) && is.finite(r)) r else 0
    }

    int_mu    <- safe_int(1)
    int_sigma <- safe_int(2)
    int_nu    <- safe_int(3)

    grad_eta_mu[i]    <- const * (int_mu    - emp_mu[i])    * d_mu[i]
    grad_eta_sigma[i] <- const * (int_sigma - emp_sigma[i]) * d_sigma[i]
    grad_eta_nu[i]    <- const * (int_nu    - emp_nu[i])    * d_nu[i]
  }

  grad_mu    <- as.vector(t(X_mu)    %*% grad_eta_mu)    / n
  grad_sigma <- as.vector(t(X_sigma) %*% grad_eta_sigma) / n
  grad_nu    <- as.vector(t(X_nu)    %*% grad_eta_nu)    / n

  grad_final <- c(grad_mu, grad_sigma, grad_nu)
  if (any(!is.finite(grad_final))) return(rep(1e10, p_total))
  return(grad_final)
}
