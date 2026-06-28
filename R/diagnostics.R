# ==============================================================================
# diagnostics.R
# Diagnostic tools for mdpde_bccg objects:
#   mdpde_diagnostics() - 4-panel quantile residual plot
#   mdpde_wp()          - robust worm plot
#   mdpde_weights()     - weights vs residuals scatter
# ==============================================================================

# ------------------------------------------------------------------------------
# Internal: add observation labels (via ggrepel) to an existing ggplot
# ------------------------------------------------------------------------------

#' Add index labels to outlying observations
#'
#' @param p A `ggplot` object.
#' @param df Data frame with columns `Index`, `Weights`, and the x/y columns.
#' @param label_arg Labelling argument: `"all"`, `"none"`, `FALSE`, or a
#'   numeric vector of indices to label.
#' @param weight_cutoff Weights below this threshold are considered outliers.
#' @param col_x,col_y Column names for the x and y aesthetics.
#' @return Modified `ggplot` object.
#' @keywords internal
add_labels <- function(p, df, label_arg, weight_cutoff, col_x, col_y) {
  if (is.logical(label_arg)) label_arg <- if (label_arg) "all" else "none"

  indices_to_label <- if (is.character(label_arg) && label_arg[1] == "all") {
    df$Index[df$Weights < weight_cutoff]
  } else if (is.numeric(label_arg)) {
    intersect(label_arg, df$Index)
  } else {
    return(p)
  }

  if (length(indices_to_label) == 0) return(p)

  df_labels <- df[df$Index %in% indices_to_label, ]
  p + ggrepel::geom_text_repel(
    data        = df_labels,
    ggplot2::aes(x = .data[[col_x]], y = .data[[col_y]], label = Index),
    size        = 2.75,
    box.padding = 0.5,
    max.overlaps = 20,
    color       = "black",
    fontface    = "bold",
    inherit.aes = FALSE
  )
}

# ------------------------------------------------------------------------------
# 1. Four-panel residual diagnostic plot
# ------------------------------------------------------------------------------

#' Four-panel Diagnostic Plot for an mdpde_bccg Fit
#'
#' Produces four complementary diagnostic graphics based on the randomised
#' quantile residuals of the fitted model:
#'
#' 1. **Residuals vs Fitted Values** – checks homogeneity of variance.
#' 2. **Residuals vs Index** – checks independence / sequential trends.
#' 3. **Kernel Density of Residuals** – checks normality of the residuals.
#' 4. **Normal Q-Q Plot** – checks normality; the reference line is fitted
#'    robustly from the inter-quartile range to avoid distortion by outliers.
#'
#' Under a correctly specified model, quantile residuals follow N(0, 1).
#' Points coloured red have robust weights below `weight_cutoff` and are
#' flagged as potential outliers by the estimation procedure.
#'
#' @param model An object of class `"mdpde_bccg"`.
#' @param label Labelling of potential outliers. One of:
#'   - `"all"` (default): label all observations with `weight < weight_cutoff`.
#'   - `"none"` or `FALSE`: no labels.
#'   - A numeric vector of observation indices to label.
#' @param weight_cutoff Weight threshold below which an observation is
#'   considered a potential outlier (default 0.05).
#'
#' @return Invisibly, a list with `residuals` and `weights`. The function is
#'   called for its side effect (a combined `gridExtra` plot).
#'
#' @seealso [mdpde_wp()], [mdpde_weights()], [mdpde_bccg()]
#'
#' @examples
#' \dontrun{
#' fit <- mdpde_bccg(y ~ x, data = dat, alpha = 0.15)
#' mdpde_diagnostics(fit)
#' }
#'
#' @export
mdpde_diagnostics <- function(model, label = "all", weight_cutoff = 0.05) {

  y       <- model$y
  res_q   <- model$residuals
  weights <- model$weights
  mu      <- model$fitted.values$mu
  n       <- length(y)

  df_diag <- data.frame(Index = seq_len(n), Fitted = mu,
                        Residuals = res_q, Weights = weights)
  is_out  <- Weights <- Residuals <- Fitted <- Index <- Theoretical <- NULL  # R CMD CHECK

  # 1. Residuals vs Fitted
  p1 <- ggplot2::ggplot(df_diag, ggplot2::aes(x = Fitted, y = Residuals)) +
    ggplot2::geom_point(ggplot2::aes(color = Weights < weight_cutoff), alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), guide = "none") +
    ggplot2::labs(title = "Residuals vs Fitted Values",
                  x = "Fitted Values", y = "Quantile Residuals") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # 2. Residuals vs Index
  p2 <- ggplot2::ggplot(df_diag, ggplot2::aes(x = Index, y = Residuals)) +
    ggplot2::geom_point(ggplot2::aes(color = Weights < weight_cutoff), alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), guide = "none") +
    ggplot2::labs(title = "Residuals vs Index",
                  x = "Observation index", y = "Quantile Residuals") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # 3. Kernel density
  dens_est <- density(res_q)
  df_dens  <- data.frame(x = dens_est$x, y = dens_est$y)
  x <- y_dens <- NULL
  p3 <- ggplot2::ggplot(df_dens, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line(color = "black", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    ggplot2::geom_rug(data = df_diag, ggplot2::aes(x = Residuals),
                      color = "red", sides = "b", alpha = 0.7,
                      linewidth = 0.5,
                      length = ggplot2::unit(0.04, "npc"),
                      inherit.aes = FALSE) +
    ggplot2::labs(title = "Kernel Density of Residuals",
                  x = "Quantile Residuals", y = "Density") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # 4. Q-Q plot (robust reference line from IQR)
  df_qq              <- df_diag[order(df_diag$Residuals), ]
  df_qq$Theoretical  <- qnorm(ppoints(n))

  p4 <- ggplot2::ggplot(df_qq, ggplot2::aes(x = Theoretical, y = Residuals)) +
    ggplot2::geom_point(ggplot2::aes(color = Weights < weight_cutoff), alpha = 0.6) +
    ggplot2::geom_line(ggplot2::aes(sample = Residuals), stat = "qq_line",
                       line.p = c(0.25, 0.75), color = "red", linewidth = 1) +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), guide = "none") +
    ggplot2::labs(title = "Normal Q-Q Plot",
                  x = "Theoretical quantiles", y = "Sample quantiles") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # Add labels
  p1 <- add_labels(p1, df_diag, label, weight_cutoff, "Fitted",      "Residuals")
  p2 <- add_labels(p2, df_diag, label, weight_cutoff, "Index",       "Residuals")
  p4 <- add_labels(p4, df_qq,   label, weight_cutoff, "Theoretical", "Residuals")

  gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
  invisible(list(residuals = res_q, weights = weights))
}

# ------------------------------------------------------------------------------
# 2. Robust worm plot
# ------------------------------------------------------------------------------

#' Robust Worm Plot for an mdpde_bccg Fit
#'
#' A detrended Q-Q plot (van Buuren & Fredriks, 2001) that displays the
#' deviation of each quantile residual from its theoretical N(0,1) counterpart.
#' Under a correctly specified model, all points should fluctuate randomly near
#' zero within the 95% pointwise confidence bands.
#'
#' The cubic smooth is fitted with observation weights equal to the robust MDPDE
#' weights, so that outliers do not distort the polynomial trend.
#'
#' @param model An object of class `"mdpde_bccg"`.
#' @param conf_level Confidence level for the pointwise envelope (default
#'   0.95).
#' @param label,weight_cutoff See [mdpde_diagnostics()].
#'
#' @return A `ggplot` object (invisibly returned after printing).
#'
#' @references
#' van Buuren, S., & Fredriks, M. (2001). Worm plot: a simple diagnostic device
#' for modelling growth reference curves.
#' *Statistics in Medicine*, **20**(8), 1259–1277.
#'
#' @seealso [mdpde_diagnostics()], [mdpde_weights()]
#'
#' @examples
#' \dontrun{
#' fit <- mdpde_bccg(y ~ x, data = dat, alpha = 0.15)
#' mdpde_wp(fit)
#' }
#'
#' @export
mdpde_wp <- function(model, conf_level = 0.95, label = "all",
                     weight_cutoff = 0.05) {
  n        <- model$n
  res_obs  <- model$residuals
  weights  <- model$weights

  df_raw  <- data.frame(Index = seq_len(n), Res = res_obs, Weights = weights)
  df_ord  <- df_raw[order(df_raw$Res), ]
  p_i     <- ((seq_len(n)) - 0.5) / n
  x_theo  <- qnorm(p_i)

  df_plot <- data.frame(
    Index       = df_ord$Index,
    Theoretical = x_theo,
    Deviation   = df_ord$Res - x_theo,
    Weights     = df_ord$Weights
  )

  # Pointwise 95% envelope
  z_crit   <- qnorm(1 - (1 - conf_level) / 2)
  x_grid   <- seq(-4, 4, length.out = 500)
  p_grid   <- pnorm(x_grid)
  se_grid  <- (1 / sqrt(n)) * (sqrt(p_grid * (1 - p_grid)) / dnorm(x_grid))
  df_env   <- data.frame(Theoretical = x_grid,
                          Upper = z_crit * se_grid,
                          Lower = -z_crit * se_grid)

  Theoretical <- Deviation <- Upper <- Lower <- NULL  # R CMD CHECK

  g <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = df_env,
      ggplot2::aes(x = Theoretical, ymin = Lower, ymax = Upper),
      fill = "gray85", alpha = 0.5
    ) +
    ggplot2::geom_line(data = df_env,
                       ggplot2::aes(x = Theoretical, y = Upper),
                       linetype = "dotted", color = "gray40") +
    ggplot2::geom_line(data = df_env,
                       ggplot2::aes(x = Theoretical, y = Lower),
                       linetype = "dotted", color = "gray40") +
    ggplot2::geom_hline(yintercept = 0, color = "red",
                        linewidth = 0.8, linetype = "dashed") +
    ggplot2::geom_point(
      data = df_plot,
      ggplot2::aes(x = Theoretical, y = Deviation,
                   color = Weights < weight_cutoff),
      alpha = 0.7, size = 2
    ) +
    ggplot2::geom_smooth(
      data    = df_plot,
      ggplot2::aes(x = Theoretical, y = Deviation, weight = Weights),
      method  = "lm",
      formula = y ~ poly(x, 3),
      se      = FALSE, color = "blue", linewidth = 0.8
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "black", "TRUE" = "red"), guide = "none"
    ) +
    ggplot2::labs(title = "Robust Worm Plot",
                  x = "Theoretical quantiles", y = "Deviation") +
    ggplot2::coord_cartesian(xlim = c(-4, 4), ylim = c(-2.5, 2.5)) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  g <- add_labels(g, df_plot, label, weight_cutoff,
                  col_x = "Theoretical", col_y = "Deviation")
  g
}

# ------------------------------------------------------------------------------
# 3. Weights vs residuals plot
# ------------------------------------------------------------------------------

#' Robust Weights vs Quantile Residuals Plot
#'
#' Plots the MDPDE observation weights \eqn{w_i = f(y_i; \hat\theta)^\alpha}
#' against the quantile residuals. Under a well-fitted model with outliers, two
#' clusters are expected:
#'
#' - **Normal observations** (high weight, small |residual|) grouped near
#'   \eqn{w_i \approx 1} and \eqn{r_i \approx 0}.
#' - **Potential outliers** (low weight, large |residual|) near the bottom of
#'   the plot.
#'
#' @param model An object of class `"mdpde_bccg"`.
#' @param label,weight_cutoff See [mdpde_diagnostics()].
#'
#' @return A `ggplot` object.
#'
#' @seealso [mdpde_diagnostics()], [mdpde_wp()]
#'
#' @examples
#' \dontrun{
#' fit <- mdpde_bccg(y ~ x, data = dat, alpha = 0.15)
#' mdpde_weights(fit)
#' }
#'
#' @export
mdpde_weights <- function(model, label = "all", weight_cutoff = 0.05) {

  res_q   <- model$residuals
  weights <- model$weights
  n       <- model$n

  df_plot <- data.frame(
    Index     = seq_len(n),
    Residuals = res_q,
    Weights   = weights,
    Type      = ifelse(weights < weight_cutoff, "Outlier", "Normal")
  )

  Residuals <- Weights <- Type <- NULL  # R CMD CHECK

  g <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data  = df_plot,
      ggplot2::aes(x = Residuals, y = Weights, color = Type),
      size  = 3, alpha = 0.7
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::scale_color_manual(
      values = c("Normal" = "darkblue", "Outlier" = "red")
    ) +
    ggplot2::labs(
      title = "Weights vs Quantile Residuals",
      x     = "Quantile Residuals",
      y     = "Weights"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position   = "bottom",
                   plot.title        = ggplot2::element_text(hjust = 0.5))

  g <- add_labels(g, df_plot, label, weight_cutoff,
                  col_x = "Residuals", col_y = "Weights")
  g
}
