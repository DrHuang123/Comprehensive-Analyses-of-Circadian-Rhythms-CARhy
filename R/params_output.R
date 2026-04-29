#' Fit the model and extract rhythm parameters
#'
#' @param expr_mat A numeric matrix of expression values with genes in rows and
#'   samples in columns.
#' @param time_vec A numeric vector of sampling times.
#' @param period The period of the rhythm. The default is 24.
#' @param conf_level Confidence level for confidence intervals. The default is
#'   0.95.
#'
#' @return
#' A data frame with one row per gene and the following columns:
#' \describe{
#'   \item{Gene}{Gene name.}
#'   \item{mesor}{Estimated mesor.}
#'   \item{mesor_ci_lower}{Lower bound of the confidence interval for mesor.}
#'   \item{mesor_ci_upper}{Upper bound of the confidence interval for mesor.}
#'   \item{amplitude}{Estimated amplitude.}
#'   \item{amplitude_ci_lower}{Lower bound of the approximate confidence
#'     interval for amplitude.}
#'   \item{amplitude_ci_upper}{Upper bound of the approximate confidence
#'     interval for amplitude.}
#'   \item{phase}{Estimated phase in the same unit as \code{time_vec}.}
#'   \item{phase_ci_lower}{Lower bound of the approximate confidence interval
#'     for phase.}
#'   \item{phase_ci_upper}{Upper bound of the approximate confidence interval
#'     for phase.}
#' }
#'
#' @details
#' For each gene, missing values are removed together with their corresponding
#' sampling times before model fitting.
#'
#' When the estimated amplitude was numerically close to zero, the phase estimate
#' was not considered well defined and the amplitude confidence interval was
#' considered numerically unstable. In this case, phase, amplitude confidence
#' interval, and phase confidence interval are returned as NA.
#'
#' Phase and phase confidence interval values may be negative. For example,
#' when \code{period = 24}, a phase of \code{-2} h is equivalent to \code{22} h.
#'
#' @references
#' Oehlert, G. W. (1992).
#' A note on the delta method. \emph{The American Statistician}, 46(1), 27--29.
#'
#' Huang, W., Menet, J., and Sinha, S. (2026).
#' CARhy: Comprehensive Analyses of Circadian Rhythms in Transcriptomic
#' Experiments with Multiple Conditions.
#'
#' @examples
#' # ------------------------------------------------------------
#' # Standard example:
#' # evenly spaced sampling times with equal numbers of
#' # replicates across time points
#' # ------------------------------------------------------------
#' set.seed(1)
#'
#' period <- 24
#' time_vec <- rep(c(2, 6, 10, 14, 18, 22), each = 3)
#' amp <- c(1.2, 1.0, 0.8, 0.6, 0.4, 0.2)
#' phi <- c(2, 4, 6, 8, 10, 12)
#' mesor <- 1
#'
#' expr_mat <- t(sapply(seq_along(amp), function(i) {
#'   mesor + amp[i] * cos(2 * pi * (time_vec - phi[i]) / period) +
#'     rnorm(length(time_vec), 0, 0.2)
#' }))
#'
#' rownames(expr_mat) <- paste0("gene", seq_along(amp))
#' colnames(expr_mat) <- paste0(
#'   "Time ", time_vec, ", rep ", rep(1:3, times = 6)
#' )
#'
#' param_res <- params_output(expr_mat, time_vec, period = 24)
#' param_res
#'
#' # ------------------------------------------------------------
#' # Nonstandard example:
#' # unevenly spaced sampling times with unequal numbers of
#' # replicates across time points and missing values
#' # ------------------------------------------------------------
#' set.seed(1)
#'
#' period <- 24
#' time_vec2 <- c(
#'   rep(2, 1),
#'   rep(6, 1),
#'   rep(10, 2),
#'   rep(12, 2),
#'   rep(18, 3),
#'   rep(22, 3)
#' )
#'
#' amp2 <- c(1.2, 1.0, 0.8, 0.6, 0.4, 0.2)
#' phi2 <- c(2, 4, 6, 8, 10, 12)
#' mesor2 <- 1
#'
#' expr_mat2 <- t(sapply(seq_along(amp2), function(i) {
#'   mesor2 + amp2[i] * cos(2 * pi * (time_vec2 - phi2[i]) / period) +
#'     rnorm(length(time_vec2), 0, 0.2)
#' }))
#'
#' rownames(expr_mat2) <- paste0("gene", seq_along(amp2))
#' colnames(expr_mat2) <- paste0(
#'   "Time ", time_vec2, ", rep ",
#'   ave(time_vec2, time_vec2, FUN = seq_along)
#' )
#'
#' # For gene5 and gene6, keep only one observation at 10 h and 12 h;
#' # set the remaining observations at those time points to NA.
#' # time 10 h -> columns 3, 4
#' # time 12 h -> columns 5, 6
#' expr_mat2[c("gene5", "gene6"), c(4, 6)] <- NA
#'
#' param_res2 <- params_output(expr_mat2, time_vec2, period = 24)
#' param_res2
#'
#' @export

params_output <- function(expr_mat, time_vec, period = 24, conf_level = 0.95) {
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "double"

  if (!is.numeric(time_vec)) {
    stop("`time_vec` must be numeric.")
  }

  if (!is.numeric(period) || length(period) != 1 ||
      !is.finite(period) || period <= 0) {
    stop("`period` must be a single positive numeric value.")
  }

  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single numeric value between 0 and 1.")
  }

  if (ncol(expr_mat) != length(time_vec)) {
    stop("`ncol(expr_mat)` must equal `length(time_vec)`.")
  }

  gene_names <- rownames(expr_mat)
  if (is.null(gene_names)) {
    gene_names <- paste0("gene", seq_len(nrow(expr_mat)))
  }

  w <- 2 * pi / period
  n_gene <- nrow(expr_mat)
  amp_tol <- sqrt(.Machine$double.eps)

  out <- data.frame(
    Gene = gene_names,
    mesor = NA_real_,
    mesor_ci_lower = NA_real_,
    mesor_ci_upper = NA_real_,
    amplitude = NA_real_,
    amplitude_ci_lower = NA_real_,
    amplitude_ci_upper = NA_real_,
    phase = NA_real_,
    phase_ci_lower = NA_real_,
    phase_ci_upper = NA_real_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  for (i in seq_len(n_gene)) {
    y_all <- as.numeric(expr_mat[i, ])
    ok <- is.finite(y_all) & is.finite(time_vec)

    y <- y_all[ok]
    t <- time_vec[ok]

    if (length(y) <= 3) {
      next
    }

    if (length(unique(t)) < 3) {
      next
    }

    X <- cbind(
      intercept = 1,
      cos = cos(w * t),
      sin = sin(w * t)
    )

    if (qr(X)$rank < 3) {
      next
    }

    XtXinv <- tryCatch(
      solve(crossprod(X)),
      error = function(e) NULL
    )

    if (is.null(XtXinv)) {
      next
    }

    beta <- XtXinv %*% crossprod(X, y)

    fitted <- as.numeric(X %*% beta)
    resid <- y - fitted
    df <- length(y) - ncol(X)

    if (df <= 0) {
      next
    }

    sigma2 <- sum(resid^2) / df

    if (!is.finite(sigma2) || sigma2 < 0) {
      next
    }

    cov_beta <- sigma2 * XtXinv

    b0 <- as.numeric(beta[1])
    b1 <- as.numeric(beta[2])
    b2 <- as.numeric(beta[3])

    mesor_hat <- b0
    amp_hat <- sqrt(b1^2 + b2^2)

    alpha <- 1 - conf_level
    tcrit <- stats::qt(1 - alpha / 2, df = df)

    mesor_ci_lower <- NA_real_
    mesor_ci_upper <- NA_real_

    amp_ci_lower <- NA_real_
    amp_ci_upper <- NA_real_

    phase_hat <- NA_real_
    phase_ci_lower <- NA_real_
    phase_ci_upper <- NA_real_

    if (is.finite(cov_beta[1, 1]) && cov_beta[1, 1] >= 0) {
      mesor_se_internal <- sqrt(cov_beta[1, 1])
      mesor_ci_lower <- mesor_hat - tcrit * mesor_se_internal
      mesor_ci_upper <- mesor_hat + tcrit * mesor_se_internal
    }

    if (is.finite(amp_hat) && amp_hat > amp_tol) {
      denom <- b1^2 + b2^2

      grad_log_amp <- c(
        0,
        b1 / denom,
        b2 / denom
      )

      log_amp_var <- as.numeric(
        t(grad_log_amp) %*% cov_beta %*% grad_log_amp
      )

      if (is.finite(log_amp_var) && log_amp_var >= 0) {
        log_amp_se_internal <- sqrt(log_amp_var)

        log_amp_lower <- log(amp_hat) - tcrit * log_amp_se_internal
        log_amp_upper <- log(amp_hat) + tcrit * log_amp_se_internal

        amp_ci_lower_tmp <- exp(log_amp_lower)
        amp_ci_upper_tmp <- exp(log_amp_upper)

        if (is.finite(amp_ci_lower_tmp) && is.finite(amp_ci_upper_tmp)) {
          amp_ci_lower <- amp_ci_lower_tmp
          amp_ci_upper <- amp_ci_upper_tmp
        }
      }

      phase_hat <- atan2(b2, b1) / w

      grad_phase <- c(
        0,
        -b2 / denom / w,
        b1 / denom / w
      )

      phase_var <- as.numeric(
        t(grad_phase) %*% cov_beta %*% grad_phase
      )

      if (is.finite(phase_var) && phase_var >= 0) {
        phase_se_internal <- sqrt(phase_var)

        phase_ci_lower <- phase_hat - tcrit * phase_se_internal
        phase_ci_upper <- phase_hat + tcrit * phase_se_internal
      }
    }

    out$mesor[i] <- mesor_hat
    out$mesor_ci_lower[i] <- mesor_ci_lower
    out$mesor_ci_upper[i] <- mesor_ci_upper

    out$amplitude[i] <- amp_hat
    out$amplitude_ci_lower[i] <- amp_ci_lower
    out$amplitude_ci_upper[i] <- amp_ci_upper

    out$phase[i] <- phase_hat
    out$phase_ci_lower[i] <- phase_ci_lower
    out$phase_ci_upper[i] <- phase_ci_upper
  }

  out
}
