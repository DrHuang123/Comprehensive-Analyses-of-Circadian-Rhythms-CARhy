#' Test the presence of rhythmicity within one condition
#'
#' @param data_list Gene-by-sample expression matrix.
#' @param time_vec Sampling-time vector.
#' @param period Expected rhythm period.
#'
#' @return A data frame with three columns: \code{pvalue} for raw p-values,
#'   \code{BH} for Benjamini-Hochberg adjusted p-values, and \code{qvalue}
#'   for Storey q-values, with one row per gene. Users may specify their own
#'   false discovery rate (FDR) level, commonly 0.05, to determine statistical
#'   significance based on \code{BH} or \code{qvalue}. The \code{qvalue} column
#'   is \code{NA} if the \pkg{qvalue} package is not installed or if q-value
#'   estimation fails.
#'
#' @references
#' Huang, W., Menet, J., and Sinha, S. (2026).
#' CARhy: Comprehensive Analyses of Circadian Rhythms in Transcriptomic Experiments with Multiple Conditions.
#'
#' @examples
#'
#' # ------------------------------------------------------------
#' # Standard example:
#' # balanced design with six evenly spaced sampling times and
#' # three replicates at each time point
#' # ------------------------------------------------------------
#' set.seed(1)
#'
#' # rhythm parameter settings
#' period <- 24
#' time_vec <- rep(c(2, 6, 10, 14, 18, 22), each = 3)
#' amp <- c(1.2, 1.0, 0.8, 0.6, 0, 0, 0.4, 0.2)          # amplitudes
#' phi <- c(2, 4, 6, 8, 0, 0, 10, 12)                    # phases
#' mesor <- 1                                            # baseline expression
#'
#' # create simulated expression data
#' data_list <- t(sapply(seq_along(amp), function(i) {
#'  mesor + amp[i] * cos(2 * pi * (time_vec - phi[i]) / period) +
#'   rnorm(length(time_vec), 0, 0.2)
#' }))
#'
#' rownames(data_list) <- paste0("gene", seq_along(amp))
#' colnames(data_list) <- paste0(
#'   "Time ", time_vec, ", rep ", rep(1:3, times = 6)
#' )
#'
#' # run TR
#' TR_res <- TR(data_list, time_vec, period)
#'
#' # view results
#' TR_res
#'
#' # count significant rhythmic genes using nominal p-values below 0.05
#' sum(TR_res$pvalue < 0.05, na.rm = TRUE)
#'
#' # specify the FDR level
#' fdr_level <- 0.05
#'
#' # count significant rhythmic genes using BH-adjusted p-values
#' sum(TR_res$BH < fdr_level, na.rm = TRUE)
#'
#' # count significant rhythmic genes using Storey q-values
#' # qvalue is NA if the qvalue package is not installed
#' sum(TR_res$qvalue < fdr_level, na.rm = TRUE)
#'
#' # ------------------------------------------------------------
#' # Nonstandard example 1:
#' # unevenly spaced sampling times (2, 6, 10, 12, 18, 22 hours)
#' # with different numbers of replicates across time points
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
#' amp2 <- c(1.2, 1.0, 0.8, 0.6, 0, 0, 0.4, 0.2)
#' phi2 <- c(2, 4, 6, 8, 0, 0, 10, 12)
#' mesor2 <- 1
#'
#' data_list2 <- t(sapply(seq_along(amp2), function(i) {
#'   mesor2 + amp2[i] * cos(2 * pi * (time_vec2 - phi2[i]) / period) +
#'     rnorm(length(time_vec2), 0, 0.2)
#' }))
#'
#' rownames(data_list2) <- paste0("gene", seq_along(amp2))
#' colnames(data_list2) <- paste0(
#'   "Time ", time_vec2, ", rep ", ave(time_vec2, time_vec2, FUN = seq_along)
#' )
#'
#' TR_res2 <- TR(data_list2, time_vec2, period)
#' TR_res2
#'
#' sum(TR_res2$pvalue < 0.05, na.rm = TRUE)
#'
#' fdr_level <- 0.05
#'
#' sum(TR_res2$BH < fdr_level, na.rm = TRUE)
#'
#' sum(TR_res2$qvalue < fdr_level, na.rm = TRUE)
#'
#' # ------------------------------------------------------------
#' # Nonstandard example 2:
#' # incomplete data with missing values at some sampling times
#' # ------------------------------------------------------------
#' set.seed(1)
#'
#' period <- 24
#' time_vec3 <- rep(c(2, 6, 10, 14, 18, 22), each = 3)
#' n_gene <- 12        # number of genes
#' amp3 <- c(1.2, 1.0, 0.8, 0.6, 0, 0, 0.4, 0.2, 1.1, 0.9, 0.7, 0.5)
#' phi3 <- c(2, 4, 6, 8, 0, 0, 10, 12, 3, 5, 7, 9)
#' mesor3 <- 1
#'
#' data_list3 <- t(sapply(seq_len(n_gene), function(i) {
#'   mesor3 + amp3[i] * cos(2 * pi * (time_vec3 - phi3[i]) / period) +
#'     rnorm(length(time_vec3), 0, 0.2)
#' }))
#'
#' rownames(data_list3) <- paste0("gene", seq_len(n_gene))
#' colnames(data_list3) <- paste0(
#'   "Time ", time_vec3, ", rep ", rep(1:3, times = 6)
#' )
#'
#' # For gene11 and gene12, keep only one observation at 2 h and 6 h;
#' # set the remaining observations at those time points to NA.
#' # time 2 h  -> columns 1, 2, 3
#' # time 6 h  -> columns 4, 5, 6
#' data_list3[c("gene11", "gene12"), c(2, 3, 5, 6)] <- NA
#'
#' TR_res3 <- TR(data_list3, time_vec3, period)
#' TR_res3
#'
#' sum(TR_res3$pvalue < 0.05, na.rm = TRUE)
#'
#' fdr_level <- 0.05
#'
#' sum(TR_res3$BH < fdr_level, na.rm = TRUE)
#'
#' sum(TR_res3$qvalue < fdr_level, na.rm = TRUE)
#'
#' @export

TR <- function(data_list, time_vec, period = 24) {

  nsimu <- 20000
  seed <- 1

  data_list <- as.matrix(data_list)

  if (ncol(data_list) != length(time_vec)) {
    stop("ncol(data_list) must equal length(time_vec).")
  }

  w <- 2 * pi / period
  G <- nrow(data_list)

  pvals <- rep(NA_real_, G)

  set.seed(seed)
  r1 <- rchisq(nsimu, 1)
  r2 <- rchisq(nsimu, 1)

  scale_cache <- list()

  for (g in seq_len(G)) {

    y_all <- as.numeric(data_list[g, ])
    ok <- is.finite(y_all) & is.finite(time_vec)

    y <- y_all[ok]
    t <- time_vec[ok]

    n <- length(y)
    df <- n - 3

    if (df <= 0) next
    if (length(unique(t)) < 3) next

    X <- cbind(1, cos(w * t), sin(w * t))

    if (qr(X)$rank < 3) next

    XtXinv <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
    if (is.null(XtXinv)) next

    V <- XtXinv[2:3, 2:3, drop = FALSE]
    eigV <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    eig1 <- eigV[1]
    eig2 <- eigV[2]

    XtY <- crossprod(X, y)
    Beta <- XtXinv %*% XtY
    Res <- as.numeric(y - X %*% Beta)
    SSE <- sum(Res^2)
    s2 <- SSE / df

    b1 <- as.numeric(Beta[2])
    b2 <- as.numeric(Beta[3])
    TS <- b1^2 + b2^2

    ev1 <- s2 * eig1
    ev2 <- s2 * eig2

    if (!(is.finite(TS) && is.finite(ev1) && is.finite(ev2) &&
          ev1 >= 0 && ev2 >= 0)) {
      next
    }

    df_key <- as.character(df)
    if (is.null(scale_cache[[df_key]])) {
      r3 <- rchisq(nsimu, df)
      scale_cache[[df_key]] <- df / r3
    }
    scale_row <- scale_cache[[df_key]]

    sim <- (r1 * ev1 + r2 * ev2) * scale_row
    pvals[g] <- mean(sim > TS)
  }

  pvals_df <- data.frame(
    pvalue = as.numeric(pvals),
    row.names = rownames(data_list)
  )

  pvals_df$BH <- stats::p.adjust(pvals_df$pvalue, method = "BH")

  pvals_df$qvalue <- NA_real_

  if (requireNamespace("qvalue", quietly = TRUE)) {
    pv <- pvals_df$pvalue
    ok_p <- is.finite(pv) & !is.na(pv) & pv >= 0 & pv <= 1
    pv_ok <- pv[ok_p]

    if (length(pv_ok) >= 2) {
      lam_max <- min(0.25, max(pv_ok) - 1e-6)

      qobj <- tryCatch(
        if (lam_max >= 0.05) {
          qvalue::qvalue(
            p = pv_ok,
            lambda = seq(0.05, lam_max, by = 0.05)
          )
        } else {
          qvalue::qvalue(
            p = pv_ok,
            lambda = 0
          )
        },
        error = function(e) {
          message("qvalue failed: ", e$message)
          NULL
        }
      )

      if (!is.null(qobj)) {
        pvals_df$qvalue[ok_p] <- qobj$qvalues
      }
    }
  }

  return(pvals_df)
}
