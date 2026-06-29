# ==============================================================================
# mdpde_bccg.R
# Main fitting function and S3 methods for class "mdpde_bccg".
# ==============================================================================

# ------------------------------------------------------------------------------
# Internal: divergence objective
# ------------------------------------------------------------------------------

#' MDPDE divergence objective for BCCG regression
#'
#' Computes the MDPDE objective function for BCCG regression:
#' \deqn{H_n(\theta) = \frac{1}{n} \sum_{i=1}^n \int f_\theta(y|x_i)^{1+\alpha} dy - \left(1+\frac{1}{\alpha}\right)\frac{1}{n}\sum_{i=1}^n f_\theta(y_i|x_i)^\alpha}
#'
#' @param theta Concatenated coefficient vector.
#' @param y,X_mu,X_sigma,X_nu Data.
#' @param alpha Tuning parameter.
#' @param mu_link,sigma_link,nu_link Link strings.
#' @return Scalar divergence value.
#' @keywords internal
#' @noRd
bccg_divergence <- function(theta, y, X_mu, X_sigma, X_nu, alpha,
                            mu_link, sigma_link, nu_link) {
  p_mu    <- ncol(X_mu); p_sigma <- ncol(X_sigma); p_nu <- ncol(X_nu)
  n       <- length(y)

  mu_i    <- link_inverse(as.vector(X_mu    %*% theta[seq_len(p_mu)]),                        mu_link)
  sigma_i <- link_inverse(as.vector(X_sigma %*% theta[(p_mu + 1):(p_mu + p_sigma)]),          sigma_link)
  nu_i    <- link_inverse(as.vector(X_nu    %*% theta[(p_mu + p_sigma + 1):(p_mu + p_sigma + p_nu)]), nu_link)

  if (any(!is.finite(mu_i)) || any(!is.finite(sigma_i)) || any(!is.finite(nu_i)) ||
      any(sigma_i <= 1e-6)  || any(mu_i <= 1e-6)) return(1e10)

  d_obs <- suppressWarnings(gamlss.dist::dBCCG(y, mu = mu_i, sigma = sigma_i, nu = nu_i))
  d_obs[!is.finite(d_obs)] <- 0
  d_obs <- pmax(d_obs, 1e-100)

  empirical_term <- (1 + 1 / alpha) * mean(d_obs^alpha)

  integrals <- numeric(n)
  for (i in seq_len(n)) {
    integrand <- function(t) {
      val <- suppressWarnings(
        gamlss.dist::dBCCG(t, mu = mu_i[i], sigma = sigma_i[i], nu = nu_i[i])
      )
      val[!is.finite(val)] <- 0
      pmax(val, 1e-100)^(1 + alpha)
    }
    res <- try(integrate(integrand, lower = 0, upper = Inf,
                         subdivisions = 100, rel.tol = 1e-2)$value, silent = TRUE)
    integrals[i] <- if (inherits(res, "try-error") || !is.finite(res)) 1e5 else res
  }

  result <- mean(integrals) - empirical_term
  if (!is.finite(result)) return(1e10)
  return(result)
}

# ------------------------------------------------------------------------------
# Internal: K matrix (empirical sandwich component)
# ------------------------------------------------------------------------------

#' Estimate the K matrix of the sandwich covariance
#'
#' K = (1/n) sum psi_i psi_i' - xi xi', where psi_i = u(y_i) f(y_i)^alpha.
#'
#' @keywords internal
calculate_k <- function(theta, y, X_mu, X_sigma, X_nu, alpha,
                        mu_link, sigma_link, nu_link) {
  p_mu    <- ncol(X_mu); p_sigma <- ncol(X_sigma); p_nu <- ncol(X_nu)
  p_total <- p_mu + p_sigma + p_nu
  n       <- length(y)

  mu_i    <- link_inverse(as.vector(X_mu    %*% theta[seq_len(p_mu)]),                        mu_link)
  sigma_i <- link_inverse(as.vector(X_sigma %*% theta[(p_mu + 1):(p_mu + p_sigma)]),          sigma_link)
  nu_i    <- link_inverse(as.vector(X_nu    %*% theta[(p_mu + p_sigma + 1):p_total]),          nu_link)

  K_matrix <- matrix(0, p_total, p_total)
  xi       <- numeric(p_total)

  for (i in seq_len(n)) {
    f_alpha_i <- gamlss.dist::dBCCG(y[i], mu = mu_i[i], sigma = sigma_i[i], nu = nu_i[i])^alpha
    u_vec     <- calculate_score(y[i], mu_i[i], sigma_i[i], nu_i[i],
                                 X_mu[i, ], X_sigma[i, ], X_nu[i, ],
                                 mu_link, sigma_link, nu_link)
    psi_i    <- u_vec * f_alpha_i
    xi       <- xi + psi_i
    K_matrix <- K_matrix + (psi_i %*% t(psi_i))
  }

  xi       <- xi / n
  K_matrix <- K_matrix / n - (xi %*% t(xi))
  return(K_matrix)
}

# ------------------------------------------------------------------------------
# Main fitting function
# ------------------------------------------------------------------------------

#' Fit a BCCG Regression Model via the MDPDE
#'
#' Estimates parameters of the Box-Cox Cole and Green (BCCG) regression model
#' using the Minimum Density Power Divergence Estimator (MDPDE) proposed by
#' Basu et al. (1998).  The tuning parameter `alpha` controls the robustness-
#' efficiency trade-off: `alpha = 0` reduces to MLE (handled via
#' [gamlss::gamlss()]), while larger values increase robustness at the cost of
#' efficiency.
#'
#' @param formula A formula for the median sub-model (e.g. `y ~ x1 + x2`).
#' @param sigma.formula A one-sided formula for the CV sub-model (default `~1`).
#' @param nu.formula A one-sided formula for the skewness sub-model
#'   (default `~1`).
#' @param alpha Non-negative tuning parameter (default 0.1).  `alpha = 0`
#'   returns an MLE fit via GAMLSS.
#' @param data A data frame.
#' @param mu.link Link function for mu: `"log"` (default), `"identity"`, or
#'   `"sqrt"`.
#' @param sigma.link Link function for sigma: `"log"` (default), `"identity"`,
#'   or `"sqrt"`.
#' @param nu.link Link function for nu: `"identity"` (default), `"log"`, or
#'   `"sqrt"`.
#'
#' @return An object of class `"mdpde_bccg"`, a list with components:
#' \describe{
#'   \item{coefficients}{Named list with entries `mu`, `sigma`, `nu`.}
#'   \item{se}{Asymptotic standard errors (sandwich), same structure.}
#'   \item{z_values, p_values}{Wald statistics and two-sided p-values.}
#'   \item{fitted.values}{Data frame with columns `mu`, `sigma`, `nu`.}
#'   \item{residuals}{Randomised quantile residuals.}
#'   \item{weights}{Robust weights \eqn{w_i = f(y_i; \hat\theta)^\alpha}.}
#'   \item{convergence}{Integer from [stats::optim()] (0 = success).}
#'   \item{cov_matrix}{Full \eqn{p \times p} sandwich covariance matrix.}
#'   \item{alpha, n, links, formula, sigma.formula, nu.formula, data, y}{
#'     Model metadata.}
#'   \item{value}{Final objective value.}
#'   \item{call}{Matched call.}
#' }
#'
#' @references
#' Basu, A., Harris, I. R., Hjort, N. L., & Jones, M. C. (1998).
#' Robust and efficient estimation by minimising a density power divergence.
#' *Biometrika*, **85**(3), 549–559. \doi{10.1093/biomet/85.3.549}
#'
#' Cole, T. J., & Green, P. J. (1992). Smoothing reference centile curves: the
#' LMS method and penalized likelihood. *Statistics in Medicine*, **11**(10),
#' 1305–1319.
#'
#' @seealso [select_alpha_mdpde()], [mdpde_diagnostics()], [summary.mdpde_bccg()]
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n <- 100
#' x <- runif(n)
#' mu_true <- exp(3 + x)
#' y <- gamlss.dist::rBCCG(n, mu = mu_true, sigma = 0.1, nu = 2)
#' dat <- data.frame(y = y, x = x)
#'
#' # Robust fit
#' fit_rob <- mdpde_bccg(y ~ x, data = dat, alpha = 0.1)
#' summary(fit_rob)
#'
#' # MLE (alpha = 0)
#' fit_mle <- mdpde_bccg(y ~ x, data = dat, alpha = 0)
#' }
#'
#' @export
mdpde_bccg <- function(formula, sigma.formula = ~1, nu.formula = ~1,
                       alpha = 0.1, data,
                       mu.link    = "log",
                       sigma.link = "log",
                       nu.link    = "identity") {
  mf      <- model.frame(formula, data = data)
  y       <- model.response(mf)
  X_mu    <- model.matrix(formula,       data = data)
  X_sigma <- model.matrix(sigma.formula, data = data)
  X_nu    <- model.matrix(nu.formula,    data = data)

  p_mu    <- ncol(X_mu); p_sigma <- ncol(X_sigma); p_nu <- ncol(X_nu)
  p_total <- p_mu + p_sigma + p_nu
  n       <- length(y)

  # ---- MLE path (alpha == 0) ------------------------------------------------
  if (alpha == 0) {
    message(">>> Alpha = 0: Returning MLE estimates via GAMLSS.")
    gamlss_call <- bquote(
      gamlss::gamlss(
        formula       = .(formula),
        sigma.formula = .(sigma.formula),
        nu.formula    = .(nu.formula),
        family        = gamlss.dist::BCCG(
          mu.link    = .(mu.link),
          sigma.link = .(sigma.link),
          nu.link    = .(nu.link)
        ),
        data  = data,
        trace = FALSE
      )
    )
    fit_mle <- try(eval(gamlss_call), silent = TRUE)
    if (inherits(fit_mle, "try-error")) stop("MLE fit (GAMLSS) failed.")

    opt <- list(
      par        = c(coef(fit_mle, "mu"), coef(fit_mle, "sigma"), coef(fit_mle, "nu")),
      convergence = as.integer(!fit_mle$converged),
      value      = -fit_mle$G.deviance / 2
    )
    cov_matrix <- vcov(fit_mle, type = "vcov")

  } else {
    # ---- Warm start ---------------------------------------------------------
    y_warm <- switch(mu.link,
                     log      = log(y),
                     sqrt     = sqrt(y),
                     identity = y)

    fit_mu_rob       <- try(robustbase::lmrob(y_warm ~ X_mu - 1), silent = TRUE)
    beta_mu_start    <- as.numeric(coef(fit_mu_rob))
    mu_hat           <- link_inverse(as.vector(X_mu %*% beta_mu_start), mu.link)

    residuals_rob    <- y - mu_hat
    sigma_cv_rob     <- mad(residuals_rob) / median(abs(mu_hat))
    sigma_start_val  <- if (sigma.link == "log") log(sigma_cv_rob) else sigma_cv_rob
    beta_sigma_start <- c(sigma_start_val, rep(0, p_sigma - 1))

    # Profile search for nu intercept
    nu_candidates    <- seq(-4, 4, by = 1)
    n_keep_nu        <- floor(n * 0.95)
    lts_errors_nu    <- vapply(nu_candidates, function(v) {
      z <- if (abs(v) < 1e-4) {
        (1 / sigma_cv_rob) * log(y / mu_hat)
      } else {
        (1 / (sigma_cv_rob * v)) * ((y / mu_hat)^v - 1)
      }
      mean(sort(z^2)[seq_len(n_keep_nu)])
    }, numeric(1))
    nu_start_val     <- nu_candidates[which.min(lts_errors_nu)]
    beta_nu_start    <- c(nu_start_val, rep(0, p_nu - 1))
    theta_start      <- c(beta_mu_start, beta_sigma_start, beta_nu_start)

    message(sprintf(">>> Alpha = %.2f: Starting robust MDPDE optimization.", alpha))

    l_sigma <- if (sigma.link == "log") log(1e-4) else 1e-4
    lower_b <- c(rep(-Inf, p_mu), rep(l_sigma, p_sigma), rep(-Inf, p_nu))

    opt <- optim(
      par     = theta_start,
      fn      = bccg_divergence,
      gr      = mdpde_gradient,
      y       = y, X_mu = X_mu, X_sigma = X_sigma, X_nu = X_nu,
      alpha   = alpha,
      mu_link = mu.link, sigma_link = sigma.link, nu_link = nu.link,
      method  = "L-BFGS-B",
      lower   = lower_b,
      control = list(
        maxit    = 1000,
        factr    = 1e7,
        pgtol    = 1e-7,
        parscale = c(rep(1, p_mu), pmax(abs(beta_sigma_start), 1), rep(2, p_nu))
      ),
      hessian = TRUE
    )

    J          <- opt$hessian
    K          <- calculate_k(opt$par, y, X_mu, X_sigma, X_nu, alpha,
                              mu.link, sigma.link, nu.link)
    cov_matrix <- tryCatch({
      J_inv <- solve(J)
      (J_inv %*% K %*% J_inv) / n
    }, error = function(e) {
      J_inv <- MASS::ginv(J)
      (J_inv %*% K %*% J_inv) / n
    })
  }

  # ---- Extract estimates ----------------------------------------------------
  theta_hat   <- opt$par
  coef_mu     <- theta_hat[seq_len(p_mu)];                         names(coef_mu)    <- colnames(X_mu)
  coef_sigma  <- theta_hat[(p_mu + 1):(p_mu + p_sigma)];           names(coef_sigma) <- colnames(X_sigma)
  coef_nu     <- theta_hat[(p_mu + p_sigma + 1):p_total];          names(coef_nu)    <- colnames(X_nu)

  se_all      <- sqrt(abs(diag(cov_matrix)))
  se_mu       <- se_all[seq_len(p_mu)];                             names(se_mu)      <- names(coef_mu)
  se_sigma    <- se_all[(p_mu + 1):(p_mu + p_sigma)];              names(se_sigma)   <- names(coef_sigma)
  se_nu       <- se_all[(p_mu + p_sigma + 1):p_total];             names(se_nu)      <- names(coef_nu)

  z_mu        <- coef_mu    / se_mu
  z_sigma     <- coef_sigma / se_sigma
  z_nu        <- coef_nu    / se_nu
  p_mu_val    <- 2 * pnorm(abs(z_mu), lower.tail = FALSE)
  p_sigma_val <- 2 * pnorm(abs(z_sigma), lower.tail = FALSE)
  p_nu_val    <- 2 * pnorm(abs(z_nu), lower.tail = FALSE)

  fitted_mu    <- link_inverse(as.vector(X_mu    %*% coef_mu),    mu.link)
  fitted_sigma <- link_inverse(as.vector(X_sigma %*% coef_sigma), sigma.link)
  fitted_nu    <- link_inverse(as.vector(X_nu    %*% coef_nu),    nu.link)

  u_val  <- pmin(pmax(
    gamlss.dist::pBCCG(y, mu = fitted_mu, sigma = fitted_sigma, nu = fitted_nu),
    1e-7), 1 - 1e-7)
  res_q  <- qnorm(u_val)

  dens_val <- gamlss.dist::dBCCG(y, mu = fitted_mu, sigma = fitted_sigma, nu = fitted_nu)
  weights  <- dens_val^alpha

  structure(
    list(
      coefficients   = list(mu = coef_mu, sigma = coef_sigma, nu = coef_nu),
      se             = list(mu = se_mu,   sigma = se_sigma,   nu = se_nu),
      z_values       = list(mu = z_mu,    sigma = z_sigma,    nu = z_nu),
      p_values       = list(mu = p_mu_val, sigma = p_sigma_val, nu = p_nu_val),
      fitted.values  = data.frame(mu = fitted_mu, sigma = fitted_sigma, nu = fitted_nu),
      residuals      = res_q,
      weights        = weights,
      convergence    = opt$convergence,
      formula        = formula,
      sigma.formula  = sigma.formula,
      nu.formula     = nu.formula,
      alpha          = alpha,
      n              = n,
      cov_matrix     = cov_matrix,
      data           = data,
      y              = y,
      links          = list(mu = mu.link, sigma = sigma.link, nu = nu.link),
      value          = opt$value,
      call           = match.call()
    ),
    class = "mdpde_bccg"
  )
}

# ------------------------------------------------------------------------------
# S3 methods
# ------------------------------------------------------------------------------

#' Print an mdpde_bccg object
#' @param x An object of class `"mdpde_bccg"`.
#' @param ... Ignored.
#' @export
print.mdpde_bccg <- function(x, ...) {
  cat("MDPDE fit for BCCG regression\nCall: ")
  print(x$call)
  cat(sprintf(
    "\nAlpha: %.2f | Divergence: %.4f | Converged: %s | n: %d\n\n",
    x$alpha, x$value,
    ifelse(x$convergence == 0, "Yes", "No"),
    x$n
  ))
  cat("Mu coefficients (link =", x$links$mu, "):\n")
  print(round(x$coefficients$mu, 4))
  cat("\nSigma coefficients (link =", x$links$sigma, "):\n")
  print(round(x$coefficients$sigma, 4))
  cat("\nNu coefficients (link =", x$links$nu, "):\n")
  print(round(x$coefficients$nu, 4))
  cat("\n(Use summary() for inference table)\n")
  invisible(x)
}

#' Summarise an mdpde_bccg object
#'
#' Prints coefficient tables with standard errors, Wald z-statistics, and p-values
#' for all three sub-models.
#'
#' @param object An object of class `"mdpde_bccg"`.
#' @param ... Ignored.
#' @export
summary.mdpde_bccg <- function(object, ...) {
  cat("MDPDE fit for BCCG regression\n")
  cat("Call: "); print(object$call)
  cat(sprintf("\nAlpha: %.2f  |  Divergence: %.4f  |  Converged: %s\n",
              object$alpha, object$value,
              ifelse(object$convergence == 0, "Yes", "No")))

  .print_table <- function(label, coef_v, se_v, z_v, p_v, link) {
    cat(sprintf(
      "\n-------------------------------------------------------------------\n%s coefficients  (link = %s)\n-------------------------------------------------------------------\n",
      toupper(label), link
    ))
    tab        <- data.frame(
      Estimate  = round(coef_v, 4),
      Std.Error = round(se_v,   7),
      z.value   = round(z_v,    4),
      p.value   = format.pval(p_v, digits = 4, eps = 1e-4)
    )
    tab$Signif <- ifelse(p_v < 0.001, "***",
                  ifelse(p_v < 0.01,  "**",
                  ifelse(p_v < 0.05,  "*",
                  ifelse(p_v < 0.1,   ".", ""))))
    print(tab)
  }

  .print_table("mu",    object$coefficients$mu,    object$se$mu,    object$z_values$mu,    object$p_values$mu,    object$links$mu)
  .print_table("sigma", object$coefficients$sigma, object$se$sigma, object$z_values$sigma, object$p_values$sigma, object$links$sigma)
  .print_table("nu",    object$coefficients$nu,    object$se$nu,    object$z_values$nu,    object$p_values$nu,    object$links$nu)
  cat("===================================================================\n")
  invisible(object)
}

#' Extract coefficients from an mdpde_bccg object
#'
#' @param object An object of class `"mdpde_bccg"`.
#' @param parameter Which sub-model: `"all"` (default), `"mu"`, `"sigma"`, or
#'   `"nu"`.
#' @param ... Ignored.
#' @export
coef.mdpde_bccg <- function(object,
                             parameter = c("all", "mu", "sigma", "nu"), ...) {
  parameter <- match.arg(parameter)
  if (parameter == "all")
    return(c(object$coefficients$mu, object$coefficients$sigma, object$coefficients$nu))
  object$coefficients[[parameter]]
}

#' Predict from an mdpde_bccg object
#'
#' @param object An `"mdpde_bccg"` fit.
#' @param newdata Optional data frame for prediction. If `NULL`, fitted values
#'   are returned.
#' @param type One of `"response"` (predicted median, default), `"parameter"`
#'   (all three parameters), or `"quantile"` (prediction interval for the
#'   response).
#' @param level Confidence / prediction level (default 0.95).
#' @param se.fit Logical; whether to include a standard error column when
#'   `type = "response"` (default `TRUE`).
#' @param ... Ignored.
#' @return A numeric vector or data frame depending on `type`.
#' @export
predict.mdpde_bccg <- function(object, newdata = NULL,
                               type  = c("response", "parameter", "quantile"),
                               level = 0.95, se.fit = TRUE, ...) {
  type <- match.arg(type)

  if (is.null(newdata)) {
    X_mu    <- model.matrix(object$formula,       data = object$data)
    X_sigma <- model.matrix(object$sigma.formula, data = object$data)
    X_nu    <- model.matrix(object$nu.formula,    data = object$data)
  } else {
    X_mu    <- model.matrix(delete.response(terms(object$formula)),       data = newdata)
    X_sigma <- model.matrix(delete.response(terms(object$sigma.formula)), data = newdata)
    X_nu    <- model.matrix(delete.response(terms(object$nu.formula)),    data = newdata)
  }

  eta_mu    <- as.vector(X_mu    %*% object$coefficients$mu)
  eta_sigma <- as.vector(X_sigma %*% object$coefficients$sigma)
  eta_nu    <- as.vector(X_nu    %*% object$coefficients$nu)

  mu_pred    <- link_inverse(eta_mu,    object$links$mu)
  sigma_pred <- link_inverse(eta_sigma, object$links$sigma)
  nu_pred    <- link_inverse(eta_nu,    object$links$nu)

  if (type == "parameter")
    return(data.frame(mu = mu_pred, sigma = sigma_pred, nu = nu_pred))

  if (type == "quantile") {
    a <- (1 - level) / 2
    return(data.frame(
      fit      = mu_pred,
      lwr_pred = gamlss.dist::qBCCG(a,     mu = mu_pred, sigma = sigma_pred, nu = nu_pred),
      upr_pred = gamlss.dist::qBCCG(1 - a, mu = mu_pred, sigma = sigma_pred, nu = nu_pred)
    ))
  }

  # type == "response"
  if (!se.fit) return(mu_pred)

  p_mu    <- length(object$coefficients$mu)
  cov_mu  <- object$cov_matrix[seq_len(p_mu), seq_len(p_mu)]
  var_eta <- rowSums((X_mu %*% cov_mu) * X_mu)
  se_eta  <- sqrt(var_eta)
  z_score <- qnorm(1 - (1 - level) / 2)

  data.frame(
    fit    = mu_pred,
    se.fit = se_eta * abs(link_derivative(mu_pred, object$links$mu)),
    lwr    = link_inverse(eta_mu - z_score * se_eta, object$links$mu),
    upr    = link_inverse(eta_mu + z_score * se_eta, object$links$mu)
  )
}
