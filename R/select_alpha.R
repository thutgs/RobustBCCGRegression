# ==============================================================================
# select_alpha.R
# Data-driven selection of the MDPDE tuning parameter alpha via the
# Trimmed Mean Negative Log-Likelihood (LVNMP) criterion.
# ==============================================================================

#' Select the Optimal MDPDE Tuning Parameter Alpha
#'
#' Implements the two-phase grid-search and local-refinement procedure described
#' in Section 3.2 of the TCC, selecting the value of `alpha` that minimises the
#' Trimmed Mean Negative Log-Likelihood (LVNMP) across the data.
#'
#' **Phase 1 – Global grid search:** fits the MDPDE for every candidate in
#' `alphas` and maps the LVNMP surface.
#'
#' **Phase 2 – Local refinement:** fits a finer grid around the Phase-1 minimum
#' (step size 0.01) when `refine = TRUE`.
#'
#' **Acceptance threshold:** the optimal alpha is forced to zero (MLE) when the
#' weighted NLL improvement over MLE is below `epsilon` (Section 3.2.3).
#'
#' @param alphas Numeric vector of candidate alpha values for Phase 1
#'   (default `seq(0, 0.35, 0.05)`).
#' @param formula,sigma.formula,nu.formula Formulas passed to [mdpde_bccg()].
#' @param data A data frame.
#' @param trim Trimming proportion for the LVNMP criterion; observations with
#'   the `trim * 100`% largest negative log-likelihood contributions are
#'   excluded (default 0.05).
#' @param epsilon Minimum relative weighted NLL improvement required to prefer
#'   the robust model over MLE (default 0.05, i.e. 5%).
#' @param only_converged If `TRUE` (default), only converged fits are eligible
#'   as optimal.
#' @param mu.link,sigma.link,nu.link Link functions passed to [mdpde_bccg()].
#' @param refine Logical; whether to run Phase 2 local refinement
#'   (default `TRUE`).
#' @param seed Numeric seed for reproducibility. Default is 123. Set to `NULL` to skip setting seed.
#'
#' @return A named list with:
#' \describe{
#'   \item{fit}{The final [mdpde_bccg()] object for the optimal alpha.}
#'   \item{optimal_alpha}{The selected alpha value.}
#'   \item{theoretical_min_alpha}{The alpha that minimised LVNMP before the
#'     threshold rule.}
#'   \item{alpha_skipped}{Logical; `TRUE` if the theoretical minimum was
#'     replaced due to non-convergence.}
#'   \item{nll_improvement}{Relative weighted NLL improvement of the robust
#'     model over MLE.}
#'   \item{table}{Data frame with columns `alpha`, `TrimError`, `Type`,
#'     `Converged`, `Category`.}
#' }
#'
#' @references
#' Hadi, A. S., & Luceno, A. (1997). Maximum trimmed likelihood estimators:
#' a unified approach, examples, and algorithms.
#' *Computational Statistics & Data Analysis*, **25**(3), 251–272.
#'
#' @seealso [mdpde_bccg()], [mdpde_diagnostics()]
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 120
#' x <- runif(n)
#' y <- gamlss.dist::rBCCG(n, mu = exp(3 + x), sigma = 0.1, nu = 2)
#' dat <- data.frame(y = y, x = x)
#'
#' sel <- select_alpha_mdpde(formula = y ~ x, data = dat)
#' cat("Optimal alpha:", sel$optimal_alpha, "\n")
#' summary(sel$fit)
#' }
#'
#' @export
select_alpha_mdpde <- function(formula,
                               sigma.formula  = ~1,
                               nu.formula     = ~1,
                               data,
                               alphas         = seq(0, 0.35, 0.05),
                               trim           = 0.05,
                               epsilon        = 0.05,
                               only_converged = TRUE,
                               mu.link        = "log",
                               sigma.link     = "log",
                               nu.link        = "identity",
                               refine         = TRUE,
                               seed = 123) {
  message(sprintf(">>> ALPHA SELECTION  (LVNMP, trim = %.0f%%)\n", trim * 100))

  if (!is.null(seed)) set.seed(seed)
  mf    <- model.frame(formula, data = data)
  y_obs <- model.response(mf)
  n     <- length(y_obs)
  n_keep <- floor(n * (1 - trim))

  # Helper: fit one alpha and compute LVNMP
  process_alpha <- function(a, type_label = "Grid") {
    tryCatch({
      invisible(capture.output({
        fit <- suppressMessages(
          mdpde_bccg(
          formula       = formula,
          sigma.formula = sigma.formula,
          nu.formula    = nu.formula,
          alpha         = a,
          data          = data,
          mu.link       = mu.link,
          sigma.link    = sigma.link,
          nu.link       = nu.link
        ))
      }))

      status_conv <- if (!is.null(fit$convergence) && fit$convergence == 0) "YES" else "NO"

      mu_est    <- fit$fitted.values$mu
      sigma_est <- fit$fitted.values$sigma
      nu_est    <- fit$fitted.values$nu

      z_i <- ifelse(
        abs(nu_est) > 1e-4,
        ((y_obs / mu_est)^nu_est - 1) / (nu_est * sigma_est),
        log(y_obs / mu_est) / sigma_est
      )
      trunc_term <- pnorm(1 / (sigma_est * abs(nu_est)))
      log_dens   <- ifelse(
        abs(nu_est) > 1e-4,
        (nu_est - 1) * log(y_obs) - nu_est * log(mu_est) -
          log(sigma_est) - 0.5 * log(2 * pi) - 0.5 * z_i^2 - log(trunc_term),
        -log(y_obs) - log(sigma_est) - 0.5 * log(2 * pi) - 0.5 * z_i^2
      )
      nll_i    <- -log_dens
      lts_val  <- mean(sort(nll_i)[seq_len(n_keep)])

      list(
        df    = data.frame(alpha = a, TrimError = lts_val, Type = type_label,
                           Converged = status_conv, stringsAsFactors = FALSE),
        model = fit
      )
    }, error = function(e) {
      list(
        df    = data.frame(alpha = a, TrimError = NA, Type = type_label,
                           Converged = "FATAL_FAIL", stringsAsFactors = FALSE),
        model = NULL
      )
    })
  }

  # Phase 1: global grid
  message(">>> Phase 1: Processing initial grid...")
  raw_results <- lapply(alphas, function(a) process_alpha(a, "Grid"))
  results     <- do.call(rbind, lapply(raw_results, `[[`, "df"))
  model_list  <- stats::setNames(lapply(raw_results, `[[`, "model"),
                                 as.character(results$alpha))

  results <- results[!is.na(results$TrimError), ]
  if (nrow(results) == 0) stop("No models converged in the initial grid.")

  res_temp       <- results[order(results$TrimError, results$alpha), ]
  best_grid_alpha <- res_temp$alpha[1]

  # Phase 2: local refinement
  if (refine && length(alphas) > 1) {
    step        <- alphas[2] - alphas[1]
    fine_alphas <- seq(max(0.01, best_grid_alpha - step),
                       min(2,    best_grid_alpha + step), 0.01)
    fine_alphas <- setdiff(round(fine_alphas, 4), round(results$alpha, 4))

    if (length(fine_alphas) > 0) {
      message(sprintf(">>> Phase 2: Refining around alpha = %.2f...", best_grid_alpha))
      fine_raw <- lapply(fine_alphas, function(a) process_alpha(a, "Refined"))
      fine_df  <- do.call(rbind, lapply(fine_raw, `[[`, "df"))
      fine_df  <- fine_df[!is.na(fine_df$TrimError), ]
      if (nrow(fine_df) > 0) {
        fine_models <- stats::setNames(lapply(fine_raw, `[[`, "model"),
                                       as.character(fine_df$alpha))
        results    <- rbind(results, fine_df)
        model_list <- c(model_list, fine_models)
      }
    }
  }

  results             <- results[order(results$TrimError, results$alpha), ]
  theoretical_min     <- results$alpha[1]
  alpha_skipped       <- FALSE

  if (only_converged) {
    conv_results <- results[results$Converged == "YES", ]
    if (nrow(conv_results) > 0) {
      optimal_alpha <- conv_results$alpha[1]
      if (optimal_alpha != theoretical_min) {
        alpha_skipped <- TRUE
        message(sprintf(
          "\n[WARNING]: Theoretical minimum (alpha=%.2f) did not converge. Selecting next best (alpha=%.2f).\n",
          theoretical_min, optimal_alpha
        ))
      }
    } else {
      optimal_alpha <- theoretical_min
    }
  } else {
    optimal_alpha <- theoretical_min
  }

  # Phase 3: acceptance threshold
  nll_improvement <- 0

  if (0 %in% results$alpha && optimal_alpha != 0) {
    fit_mle <- model_list[["0"]]
    fit_rob <- model_list[[as.character(optimal_alpha)]]

    if (!is.null(fit_mle) && !is.null(fit_rob)) {
      mu_rob    <- fit_rob$fitted.values$mu
      sigma_rob <- fit_rob$fitted.values$sigma
      nu_rob    <- fit_rob$fitted.values$nu
      mu_mle    <- fit_mle$fitted.values$mu
      sigma_mle <- fit_mle$fitted.values$sigma
      nu_mle    <- fit_mle$fitted.values$nu

      dens_rob <- pmax(gamlss.dist::dBCCG(y_obs, mu = mu_rob, sigma = sigma_rob, nu = nu_rob), 1e-300)
      dens_mle <- pmax(gamlss.dist::dBCCG(y_obs, mu = mu_mle, sigma = sigma_mle, nu = nu_mle), 1e-300)

      w_i     <- dens_rob^optimal_alpha
      sum_w   <- sum(w_i)
      wnll_rob <- sum(w_i * (-log(dens_rob))) / sum_w
      wnll_mle <- sum(w_i * (-log(dens_mle))) / sum_w

      nll_improvement <- (wnll_mle - wnll_rob) / abs(wnll_mle)

      if (nll_improvement > epsilon) {
        message(sprintf(
          "\n>>> ROBUST SELECTED: weighted NLL improvement = %.2f%% > %.0f%% threshold.\n",
          nll_improvement * 100, epsilon * 100
        ))
      } else {
        message(sprintf(
          "\n>>> MLE SELECTED: weighted NLL improvement = %.2f%% <= %.0f%% threshold.\n",
          nll_improvement * 100, epsilon * 100
        ))
        optimal_alpha <- 0
      }
    }
  }

  # Final fit
  final_fit <- model_list[[as.character(optimal_alpha)]]
  if (is.null(final_fit)) {
    invisible(capture.output({
      final_fit <- mdpde_bccg(
        formula = formula, sigma.formula = sigma.formula, nu.formula = nu.formula,
        alpha = optimal_alpha, data = data,
        mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link
      )
    }))
  }

  message(sprintf("\n>>> FINAL OPTIMAL ALPHA: %.2f <<<\n", optimal_alpha))

  summary(final_fit)

  results$Category <- results$Type
  if (nrow(results) > 0)
    results$Category[results$alpha == optimal_alpha] <- "Optimal"
  results <- results[order(results$alpha), ]

  list(
    fit                  = final_fit,
    optimal_alpha        = optimal_alpha,
    theoretical_min_alpha = theoretical_min,
    alpha_skipped        = alpha_skipped,
    nll_improvement      = nll_improvement,
    table                = results[, c("alpha", "TrimError", "Type", "Converged", "Category")]
  )
}
