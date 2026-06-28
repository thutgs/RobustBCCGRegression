# RobustBCCGRegression

Robust regression for the **Box-Cox Cole and Green (BCCG)** distribution via
the **Minimum Density Power Divergence Estimator (MDPDE)**.

This package is the computational implementation accompanying the undergraduate
thesis *"Modelos de Regressão Box-Cox Cole e Green Robustos"* (UnB, 2026).

---

## Installation

```r
# Install from GitHub:
devtools::install_github("thutgs/RobustBCCGRegression")
```

---

## Quick start

```r
library(RobustBCCGRegression)

set.seed(42)
n  <- 150
x  <- runif(n)
y  <- gamlss.dist::rBCCG(n, mu = exp(3 + x), sigma = 0.1, nu = 2)

# Introduce 5% contamination
idx_up              <- order(x)[1:4]
idx_dn              <- order(x, decreasing = TRUE)[1:4]
y[idx_up]           <- 120
y[idx_dn]           <- 0.001
dat                 <- data.frame(y = y, x = x)

# Automatic alpha selection
sel <- select_alpha_mdpde(y ~ x, data = dat)
cat("Optimal alpha:", sel$optimal_alpha, "\n")
summary(sel$fit)

# Diagnostic plots
mdpde_diagnostics(sel$fit)
mdpde_wp(sel$fit)
mdpde_weights(sel$fit)
```

---

## Key functions

| Function | Description |
|---|---|
| `mdpde_bccg()` | Fit BCCG regression for a **fixed** alpha |
| `select_alpha_mdpde()` | **Automatic** data-driven alpha selection (LVNMP) |
| `mdpde_diagnostics()` | 4-panel quantile residual diagnostic plot |
| `mdpde_wp()` | Robust worm plot |
| `mdpde_weights()` | Weights vs residuals scatter |

---

## The MDPDE in brief

The MDPDE (Basu et al., 1998) minimises the density power divergence

$$H_n(\theta) = \frac{1}{n}\sum_{i=1}^n \int_0^\infty f_\theta(y|x_i)^{1+\alpha}\,dy \;-\; \left(1+\frac{1}{\alpha}\right) \frac{1}{n}\sum_{i=1}^n f_\theta(y_i|x_i)^\alpha$$

indexed by a tuning parameter $\alpha \geq 0$.  
Setting $\alpha = 0$ recovers the **MLE**; larger values increase robustness
by downweighting observations with low fitted density (potential outliers).

---

## References

- Basu, A., Harris, I. R., Hjort, N. L., & Jones, M. C. (1998). Robust and
  efficient estimation by minimising a density power divergence. *Biometrika*,
  **85**(3), 549–559.
- Cole, T. J., & Green, P. J. (1992). Smoothing reference centile curves: the
  LMS method and penalized likelihood. *Statistics in Medicine*, **11**(10),
  1305–1319.
- Rigby, R. A., & Stasinopoulos, D. M. (2005). Generalized additive models for
  location, scale and shape. *JRSS-C*, **54**(3), 507–554.
