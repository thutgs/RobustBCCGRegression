#' RobustBCCGRegression: Robust BCCG Regression via MDPDE
#'
#' @description
#' Implements robust estimation of Box-Cox Cole and Green (BCCG) regression
#' models via the Minimum Density Power Divergence Estimator (MDPDE).
#'
#' The BCCG distribution models strictly positive, potentially skewed responses
#' through three parameters:
#' \itemize{
#'   \item \eqn{\mu_i > 0}: median sub-model.
#'   \item \eqn{\sigma_i > 0}: approximate coefficient of variation sub-model.
#'   \item \eqn{\nu_i \in \mathbb{R}}: skewness sub-model.
#' }
#'
#' The MDPDE (Basu et al., 1998) down-weights outliers via a tuning parameter
#' \eqn{\alpha \geq 0}: when \eqn{\alpha = 0} the estimator coincides with
#' the MLE (handled via \pkg{gamlss}); larger values increase robustness at
#' the cost of efficiency.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{mdpde_bccg}}}{Fit a BCCG regression model for a fixed alpha.}
#'   \item{\code{\link{select_alpha_mdpde}}}{Automatic data-driven alpha selection via
#'     the trimmed negative log-likelihood (LVNMP) criterion.}
#'   \item{\code{\link{mdpde_diagnostics}}}{Four-panel quantile residual diagnostic plot.}
#'   \item{\code{\link{mdpde_wp}}}{Robust worm plot.}
#'   \item{\code{\link{mdpde_weights}}}{Weights vs quantile residuals scatter plot.}
#' }
#'
#' @references
#' Basu, A., Harris, I. R., Hjort, N. L., & Jones, M. C. (1998).
#' Robust and efficient estimation by minimising a density power divergence.
#' \emph{Biometrika}, \strong{85}(3), 549–559.
#'
#' Cole, T. J., & Green, P. J. (1992). Smoothing reference centile curves: the
#' LMS method and penalized likelihood. \emph{Statistics in Medicine},
#' \strong{11}(10), 1305–1319.
#'
#' Rigby, R. A., & Stasinopoulos, D. M. (2005). Generalized additive models for
#' location, scale and shape. \emph{Journal of the Royal Statistical Society
#' Series C}, \strong{54}(3), 507–554.
#'
#' @import gamlss
#' @import gamlss.dist
#' @import robustbase
#' @import MASS
#' @import ggplot2
#' @import ggrepel
#' @import gridExtra
#' @importFrom stats model.frame model.matrix model.response optim qnorm pnorm dnorm pt integrate coef vcov mad median density ppoints delete.response terms setNames
#' @importFrom utils capture.output
#'
#' @name RobustBCCGRegression-package
#' @aliases RobustBCCGRegression
#'
"_PACKAGE"

utils::globalVariables(c(
  "Index", "Weights", "Fitted", "Residuals",
  "Theoretical", "Upper", "Lower", "Type",
  "x", "y", "y_dens"
))
