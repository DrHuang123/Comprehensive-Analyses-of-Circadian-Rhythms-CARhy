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
#'   \item{mesor_se}{Standard error of the estimated mesor.}
#'   \item{mesor_ci_lower}{Lower bound of the confidence interval for mesor.}
#'   \item{mesor_ci_upper}{Upper bound of the confidence interval for mesor.}
#'   \item{amplitude}{Estimated amplitude.}
#'   \item{amplitude_se}{Approximate standard error of the estimated amplitude,
#'     calculated using the delta method.}
#'   \item{amplitude_ci_lower}{Lower bound of the approximate confidence
#'     interval for amplitude, truncated at 0.}
#'   \item{amplitude_ci_upper}{Upper bound of the approximate confidence
#'     interval for amplitude.}
#'   \item{phase}{Estimated phase in the same unit as \code{time_vec}.}
#'   \item{phase_se}{Approximate standard error of the estimated phase,
#'     calculated using the delta method.}
#'   \item{phase_ci_lower}{Lower bound of the approximate confidence interval
#'     for phase, wrapped to the interval \code{[0, period)}.}
#'   \item{phase_ci_upper}{Upper bound of the approximate confidence interval
#'     for phase, wrapped to the interval \code{[0, period)}.}
#' }
#'
#' @details
#' For each gene, missing values are removed together with their corresponding
#' sampling times before model fitting. The model is fitted as
#' \code{expression ~ cos(2 * pi * time / period) + sin(2 * pi * time / period)}.
#'
#' Standard errors for mesor are obtained from the fitted linear model. Standard
#' errors for amplitude and phase are approximate and are calculated using the
#' delta method. Because amplitude is restricted to be non-negative, the lower
#' confidence bound for amplitude is truncated at 0. Because phase is circular,
#' phase confidence intervals are wrapped to \code{[0, period)}. If
#' \code{phase_ci_lower} is greater than \code{phase_ci_upper}, the phase
#' confidence interval crosses the circular boundary at 0 or \code{period}.
#'
#' When the estimated amplitude is close to zero, phase is not well defined. In
#' this case, the phase standard error and confidence interval are returned as
#' \code{NA}.
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

  out <- data.frame(
    Gene = gene_names,
    mesor = NA_real_,
    mesor_se = NA_real_,
    mesor_ci_lower = NA_real_,
    mesor_ci_upper = NA_real_,
    amplitude = NA_real_,
    amplitude_se = NA_real_,
    amplitude_ci_lower = NA_real_,
    amplitude_ci_upper = NA_real_,
    phase = NA_real_,
    phase_se = NA_real_,
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

    XtXinv <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
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

    phase_hat <- atan2(b2, b1) / w
    phase_hat <- phase_hat %% period

    alpha <- 1 - conf_level
    tcrit <- stats::qt(1 - alpha / 2, df = df)

    mesor_se <- NA_real_
    mesor_ci_lower <- NA_real_
    mesor_ci_upper <- NA_real_

    if (is.finite(cov_beta[1, 1]) && cov_beta[1, 1] >= 0) {
      mesor_se <- sqrt(cov_beta[1, 1])
      mesor_ci_lower <- mesor_hat - tcrit * mesor_se
      mesor_ci_upper <- mesor_hat + tcrit * mesor_se
    }

    amp_se <- NA_real_
    amp_ci_lower <- NA_real_
    amp_ci_upper <- NA_real_

    phase_se <- NA_real_
    phase_ci_lower <- NA_real_
    phase_ci_upper <- NA_real_

    if (is.finite(amp_hat) && amp_hat > sqrt(.Machine$double.eps)) {
      grad_amp <- c(0, b1 / amp_hat, b2 / amp_hat)
      amp_var <- as.numeric(t(grad_amp) %*% cov_beta %*% grad_amp)

      if (is.finite(amp_var) && amp_var >= 0) {
        amp_se <- sqrt(amp_var)
        amp_ci_lower <- max(0, amp_hat - tcrit * amp_se)
        amp_ci_upper <- amp_hat + tcrit * amp_se
      }

      denom <- b1^2 + b2^2
      grad_phase <- c(0, -b2 / denom / w, b1 / denom / w)
      phase_var <- as.numeric(t(grad_phase) %*% cov_beta %*% grad_phase)

      if (is.finite(phase_var) && phase_var >= 0) {
        phase_se <- sqrt(phase_var)

        phase_lower_unwrapped <- phase_hat - tcrit * phase_se
        phase_upper_unwrapped <- phase_hat + tcrit * phase_se

        phase_ci_lower <- phase_lower_unwrapped %% period
        phase_ci_upper <- phase_upper_unwrapped %% period
      }
    }

    out$mesor[i] <- mesor_hat
    out$mesor_se[i] <- mesor_se
    out$mesor_ci_lower[i] <- mesor_ci_lower
    out$mesor_ci_upper[i] <- mesor_ci_upper

    out$amplitude[i] <- amp_hat
    out$amplitude_se[i] <- amp_se
    out$amplitude_ci_lower[i] <- amp_ci_lower
    out$amplitude_ci_upper[i] <- amp_ci_upper

    out$phase[i] <- phase_hat
    out$phase_se[i] <- phase_se
    out$phase_ci_lower[i] <- phase_ci_lower
    out$phase_ci_upper[i] <- phase_ci_upper
  }

  out
}
