#' Test whether mesor differs across multiple conditions.
#'
#' @param data_list List of condition-specific expression matrices.
#' @param timepoint_list List of condition-specific sampling-time vectors.
#' @param ncond Number of conditions.
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
#' time_vec1 <- rep(c(2, 6, 10, 14, 18, 22), each = 3)
#' time_vec2 <- rep(c(2, 6, 10, 14, 18, 22), each = 3)
#' time_vec3 <- rep(c(2, 6, 10, 14, 18, 22), each = 3)
#'
#' amp <- c(1.2, 1.0, 0.8, 0.6, 0.4, 0.4, 0.2, 0.2)      # amplitudes
#' phi <- c(2, 4, 6, 8, 10, 12, 14, 16)                  # phases
#'
#' # mesors for multiple conditions
#' mesor1 <- rep(1, 8)
#' mesor2 <- c(1.8, 1.6, 1.4, 1.2, 1, 1, 1, 1)
#' mesor3 <- c(2.0, 1.7, 1.5, 1.3, 1, 1, 1, 1)
#'
#'
#' # create simulated expression data
#' data1 <- t(sapply(seq_along(amp), function(i) {
#'   mesor1[i] + amp[i] * cos(2 * pi * (time_vec1 - phi[i]) / period) +
#'     rnorm(length(time_vec1), 0, 0.2)
#' }))
#'
#' data2 <- t(sapply(seq_along(amp), function(i) {
#'   mesor2[i] + amp[i] * cos(2 * pi * (time_vec2 - phi[i]) / period) +
#'     rnorm(length(time_vec2), 0, 0.2)
#' }))
#'
#' data3 <- t(sapply(seq_along(amp), function(i) {
#'   mesor3[i] + amp[i] * cos(2 * pi * (time_vec3 - phi[i]) / period) +
#'     rnorm(length(time_vec3), 0, 0.2)
#' }))
#'
#' rownames(data1) <- paste0("gene", seq_along(amp))
#' rownames(data2) <- paste0("gene", seq_along(amp))
#' rownames(data3) <- paste0("gene", seq_along(amp))
#' colnames(data1) <- paste0("Time ", time_vec1, ", rep ", rep(1:3, times = 6))
#' colnames(data2) <- paste0("Time ", time_vec2, ", rep ", rep(1:3, times = 6))
#' colnames(data3) <- paste0("Time ", time_vec3, ", rep ", rep(1:3, times = 6))
#'
#'
#' # combine data and time points as lists
#' data_list <- list(data1, data2, data3)
#' timepoint_list <- list(time_vec1, time_vec2, time_vec3)
#'
#' # run TDM
#' TDM_res <- TDM(data_list = data_list, timepoint_list = timepoint_list, ncond = 3, period = period)
#'
#' # view results
#' TDM_res
#'
#' # count genes with significant differential mesor using nominal p-values below 0.05
#' sum(TDM_res$pvalue < 0.05, na.rm = TRUE)
#'
#' # specify the FDR level
#' fdr_level <- 0.05
#'
#' # count genes with significant differential mesor using BH-adjusted p-values
#' sum(TDM_res$BH < fdr_level, na.rm = TRUE)
#'
#' # count genes with significant differential mesor using Storey q-values
#' # qvalue is NA if the qvalue package is not installed
#' sum(TDM_res$qvalue < fdr_level, na.rm = TRUE)
#'
#' # ------------------------------------------------------------
#' # Nonstandard example 1:
#' # unevenly spaced sampling times (2, 6, 10, 12, 18, 22 hours)
#' # with different numbers of replicates across time points
#' # and experimental conditions
#' # ------------------------------------------------------------
#' set.seed(1)
#'
#' period <- 24
#'
#' time_vec1_ns <- c(
#'   rep(2, 1), rep(6, 1), rep(10, 2),
#'   rep(12, 2), rep(18, 3), rep(22, 3)
#' )
#'
#' time_vec2_ns <- c(
#'   rep(2, 2), rep(6, 2), rep(10, 3),
#'   rep(12, 3), rep(18, 4), rep(22, 4)
#' )
#'
#' time_vec3_ns <- c(
#'   rep(2, 2), rep(6, 2), rep(10, 3),
#'   rep(12, 3), rep(18, 3), rep(22, 4)
#' )
#'
#' amp_ns <- c(
#'   1.2, 1.0, 0.8, 0.6, 0.4, 0.4, 0.2, 0.2, 0.9, 0.7,
#'   0.5, 0.3, 1.2, 1.0, 0.8, 0.6, 0.4, 0.4, 0.2, 0.2
#' )
#'
#' phi_ns <- c(
#'   2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
#'   22, 24, 2, 4, 6, 8, 10, 12, 14, 16
#' )
#'
#' mesor1_ns <- c(
#'   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
#'   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
#' )
#'
#' mesor2_ns <- c(
#'   1.8, 1.6, 1.4, 1.2, 1.0, 1.0, 1.0, 1.0, 1.7, 1.5,
#'   1.3, 1.1, 1.8, 1.6, 1.4, 1.2, 1.0, 1.0, 1.0, 1.0
#' )
#'
#' mesor3_ns <- c(
#'   2.0, 1.7, 1.5, 1.3, 1.0, 1.0, 1.0, 1.0, 1.9, 1.6,
#'   1.4, 1.2, 2.0, 1.7, 1.5, 1.3, 1.0, 1.0, 1.0, 1.0
#' )
#'
#' data1_ns <- t(sapply(seq_along(amp_ns), function(i) {
#'   mesor1_ns[i] + amp_ns[i] * cos(2 * pi * (time_vec1_ns - phi_ns[i]) / period) +
#'     rnorm(length(time_vec1_ns), 0, 0.2)
#' }))
#'
#' data2_ns <- t(sapply(seq_along(amp_ns), function(i) {
#'   mesor2_ns[i] + amp_ns[i] * cos(2 * pi * (time_vec2_ns - phi_ns[i]) / period) +
#'     rnorm(length(time_vec2_ns), 0, 0.2)
#' }))
#'
#' data3_ns <- t(sapply(seq_along(amp_ns), function(i) {
#'   mesor3_ns[i] + amp_ns[i] * cos(2 * pi * (time_vec3_ns - phi_ns[i]) / period) +
#'     rnorm(length(time_vec3_ns), 0, 0.2)
#' }))
#'
#' rownames(data1_ns) <- paste0("gene", seq_along(amp_ns))
#' rownames(data2_ns) <- paste0("gene", seq_along(amp_ns))
#' rownames(data3_ns) <- paste0("gene", seq_along(amp_ns))
#' colnames(data1_ns) <- paste0(
#'   "Cond1_Time ", time_vec1_ns, ", rep ",
#'   ave(time_vec1_ns, time_vec1_ns, FUN = seq_along)
#' )
#' colnames(data2_ns) <- paste0(
#'   "Cond2_Time ", time_vec2_ns, ", rep ",
#'   ave(time_vec2_ns, time_vec2_ns, FUN = seq_along)
#' )
#' colnames(data3_ns) <- paste0(
#'   "Cond3_Time ", time_vec3_ns, ", rep ",
#'   ave(time_vec3_ns, time_vec3_ns, FUN = seq_along)
#' )
#'
#' data_list_ns <- list(data1_ns, data2_ns, data3_ns)
#' timepoint_list_ns <- list(time_vec1_ns, time_vec2_ns, time_vec3_ns)
#'
#' TDM_res_ns <- TDM(
#'   data_list = data_list_ns,
#'   timepoint_list = timepoint_list_ns,
#'   ncond = 3,
#'   period = 24
#' )
#'
#' TDM_res_ns
#'
#' sum(TDM_res_ns$pvalue < 0.05, na.rm = TRUE)
#'
#' fdr_level <- 0.05
#'
#' sum(TDM_res_ns$BH < fdr_level, na.rm = TRUE)
#'
#' sum(TDM_res_ns$qvalue < fdr_level, na.rm = TRUE)
#'
#' # ------------------------------------------------------------
#' # Nonstandard example 2:
#' # incomplete data with missing values at some sampling times
#' # ------------------------------------------------------------
#' set.seed(1)
#'
#' period <- 24
#' time_vec1_miss <- rep(c(2, 6, 10, 14, 18, 22), each = 3)
#' time_vec2_miss <- rep(c(2, 6, 10, 14, 18, 22), each = 3)
#' time_vec3_miss <- rep(c(2, 6, 10, 14, 18, 22), each = 3)
#'
#' amp_miss <- c(
#'   1.2, 1.0, 0.8, 0.6, 0.4, 0.4, 0.2, 0.2, 0.9, 0.7,
#'   0.5, 0.3, 1.2, 1.0, 0.8, 0.6, 0.4, 0.4, 0.2, 0.2
#' )
#'
#' phi_miss <- c(
#'   2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
#'   22, 24, 2, 4, 6, 8, 10, 12, 14, 16
#' )
#'
#' mesor1_miss <- c(
#'   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
#'   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
#' )
#'
#' mesor2_miss <- c(
#'   1.8, 1.6, 1.4, 1.2, 1.0, 1.0, 1.0, 1.0, 1.7, 1.5,
#'   1.3, 1.1, 1.8, 1.6, 1.4, 1.2, 1.0, 1.0, 1.0, 1.0
#' )
#'
#' mesor3_miss <- c(
#'   2.0, 1.7, 1.5, 1.3, 1.0, 1.0, 1.0, 1.0, 1.9, 1.6,
#'   1.4, 1.2, 2.0, 1.7, 1.5, 1.3, 1.0, 1.0, 1.0, 1.0
#' )
#'
#' data1_miss <- t(sapply(seq_along(amp_miss), function(i) {
#'   mesor1_miss[i] + amp_miss[i] * cos(2 * pi * (time_vec1_miss - phi_miss[i]) / period) +
#'     rnorm(length(time_vec1_miss), 0, 0.2)
#' }))
#'
#' data2_miss <- t(sapply(seq_along(amp_miss), function(i) {
#'   mesor2_miss[i] + amp_miss[i] * cos(2 * pi * (time_vec2_miss - phi_miss[i]) / period) +
#'     rnorm(length(time_vec2_miss), 0, 0.2)
#' }))
#'
#' data3_miss <- t(sapply(seq_along(amp_miss), function(i) {
#'   mesor3_miss[i] + amp_miss[i] * cos(2 * pi * (time_vec3_miss - phi_miss[i]) / period) +
#'     rnorm(length(time_vec3_miss), 0, 0.2)
#' }))
#'
#' rownames(data1_miss) <- paste0("gene", seq_along(amp_miss))
#' rownames(data2_miss) <- paste0("gene", seq_along(amp_miss))
#' rownames(data3_miss) <- paste0("gene", seq_along(amp_miss))
#' colnames(data1_miss) <- paste0(
#'   "Cond1_Time ", time_vec1_miss, ", rep ",
#'   ave(time_vec1_miss, time_vec1_miss, FUN = seq_along)
#' )
#' colnames(data2_miss) <- paste0(
#'   "Cond2_Time ", time_vec2_miss, ", rep ",
#'   ave(time_vec2_miss, time_vec2_miss, FUN = seq_along)
#' )
#' colnames(data3_miss) <- paste0(
#'   "Cond3_Time ", time_vec3_miss, ", rep ",
#'   ave(time_vec3_miss, time_vec3_miss, FUN = seq_along)
#' )
#'
#' # For gene19 and gene20 in condition 1, keep only one observation
#' # at 2 h and 6 h; set the remaining observations at those time points to NA.
#' # time 2 h -> columns 1, 2, 3
#' # time 6 h -> columns 4, 5, 6
#' data1_miss[c("gene19", "gene20"), c(2, 3, 5, 6)] <- NA
#'
#' data_list_miss <- list(data1_miss, data2_miss, data3_miss)
#' timepoint_list_miss <- list(time_vec1_miss, time_vec2_miss, time_vec3_miss)
#'
#' TDM_res_miss <- TDM(
#'   data_list = data_list_miss,
#'   timepoint_list = timepoint_list_miss,
#'   ncond = 3,
#'   period = 24
#' )
#'
#' TDM_res_miss
#'
#' sum(TDM_res_miss$pvalue < 0.05, na.rm = TRUE)
#'
#' fdr_level <- 0.05
#'
#' sum(TDM_res_miss$BH < fdr_level, na.rm = TRUE)
#'
#' sum(TDM_res_miss$qvalue < fdr_level, na.rm = TRUE)
#'
#' @export

TDM <- function(data_list, timepoint_list, ncond, period = 24) {

  ncond <- length(data_list)

  if (!is.list(data_list) || length(data_list) < 2) {
    stop("`data_list` must be a list with at least two expression matrices.")
  }

  if (length(data_list) != ncond) {
    stop("`ncond` must equal `length(data_list)`.")
  }

  if (!is.list(timepoint_list)) {
    shared_time <- as.numeric(timepoint_list)
    timepoint_list <- replicate(ncond, shared_time, simplify = FALSE)
  }

  if (length(timepoint_list) != ncond) {
    stop("`timepoint_list` must have the same length as `data_list`.")
  }

  if (!is.numeric(period) || length(period) != 1L || is.na(period) || period <= 0) {
    stop("`period` must be one positive number.")
  }

  data_list <- lapply(data_list, function(x) {
    x <- as.matrix(x)
    x
  })

  n_genes <- nrow(data_list[[1]])
  if (any(vapply(data_list, nrow, integer(1)) != n_genes)) {
    stop("All matrices in `data_list` must have the same number of rows.")
  }

  timepoint_list <- lapply(timepoint_list, as.numeric)

  for (k in seq_len(ncond)) {
    if (ncol(data_list[[k]]) != length(timepoint_list[[k]])) {
      stop(
        "For condition ", k,
        ", `ncol(data_list[[k]])` must equal `length(timepoint_list[[k]])`."
      )
    }
  }

  w <- 2 * pi / period

  X_list <- vector("list", ncond)
  xtxinv_list <- vector("list", ncond)
  xtxinvxt_list <- vector("list", ncond)
  df_vec <- integer(ncond)

  for (k in seq_len(ncond)) {
    tt <- timepoint_list[[k]]
    Xk <- cbind(1, cos(w * tt), sin(w * tt))
    df_k <- nrow(Xk) - 3L

    if (df_k <= 0) {
      stop("Condition ", k, " needs at least 4 samples/time points. Got n = ", nrow(Xk))
    }

    xtx <- crossprod(Xk)
    xtxinv <- solve(xtx)

    X_list[[k]] <- Xk
    xtxinv_list[[k]] <- xtxinv
    xtxinvxt_list[[k]] <- xtxinv %*% t(Xk)
    df_vec[k] <- df_k
  }

  rho <- ncond - 1L
  Lmat <- matrix(0, nrow = rho, ncol = 3L * ncond)
  for (j in seq_len(rho)) {
    Lmat[j, 1] <- 1
    Lmat[j, 3L * j + 1L] <- -1
  }

  block_mats <- vector("list", ncond)
  for (k in seq_len(ncond)) {
    M <- matrix(0, nrow = 3L * ncond, ncol = 3L * ncond)
    idx <- (3L * (k - 1L) + 1L):(3L * k)
    M[idx, idx] <- xtxinv_list[[k]]
    block_mats[[k]] <- M
  }

  pvals <- rep(NA_real_, n_genes)

  for (g in seq_len(n_genes)) {
    K <- length(data_list)

    gamma_vec <- numeric(3L * K)
    sigma2_hat <- numeric(K)
    varcov <- matrix(0, nrow = 3L * K, ncol = 3L * K)

    bad_gene <- FALSE

    for (k in seq_len(K)) {
      y <- as.numeric(data_list[[k]][g, ])

      if (anyNA(y)) {
        bad_gene <- TRUE
        break
      }

      gamma_hat <- xtxinvxt_list[[k]] %*% y
      resid <- y - as.vector(X_list[[k]] %*% gamma_hat)
      sigma2_k <- sum(resid^2) / df_vec[k]

      if (!is.finite(sigma2_k) || sigma2_k < 0) {
        bad_gene <- TRUE
        break
      }

      idx <- (3L * (k - 1L) + 1L):(3L * k)
      gamma_vec[idx] <- as.numeric(gamma_hat)
      sigma2_hat[k] <- sigma2_k
      varcov[idx, idx] <- sigma2_k * xtxinv_list[[k]]
    }

    if (bad_gene) {
      pvals[g] <- NA_real_
      next
    }

    middle <- Lmat %*% varcov %*% t(Lmat)
    Lcov <- tryCatch(solve(middle), error = function(e) NULL)
    if (is.null(Lcov)) {
      pvals[g] <- NA_real_
      next
    }

    diff_vec <- Lmat %*% gamma_vec
    stat <- as.numeric(t(diff_vec) %*% Lcov %*% diff_vec)
    if (!is.finite(stat) || stat < 0) {
      pvals[g] <- NA_real_
      next
    }

    nrlma <- nrow(Lmat)
    mu1 <- nrlma
    mu2 <- 2 * nrlma

    store_var_sigma2 <- numeric(K)
    store_capb <- vector("list", K)

    for (k in seq_len(K)) {
      capbk <- Lmat %*% block_mats[[k]] %*% t(Lmat) %*% Lcov
      store_capb[[k]] <- capbk

      term2 <- capbk %*% capbk
      term3 <- term2 %*% capbk
      term4 <- term3 %*% capbk

      term1_trace <- sum(diag(capbk))
      term2_trace <- sum(diag(term2))
      term3_trace <- sum(diag(term3))
      term4_trace <- sum(diag(term4))

      df_k <- df_vec[k]
      var_sigma2_k <- 2 * sigma2_hat[k]^2 / df_k
      store_var_sigma2[k] <- var_sigma2_k

      mu1 <- mu1 + 2 * term2_trace * sigma2_hat[k]^2 / df_k

      mu1_raw <- df_k
      mu2_raw <- 2 * (df_k / 2 + 1) * df_k
      mu3_raw <- 4 * (df_k / 2 + 2) * (df_k / 2 + 1) * df_k
      mu4_raw <- 8 * (df_k / 2 + 3) * (df_k / 2 + 2) * (df_k / 2 + 1) * df_k

      third_cnt_sigma2_k <- (sigma2_hat[k]^3 / df_k^3) *
        (mu3_raw - 3 * mu2_raw * mu1_raw + 2 * mu1_raw^3)

      fourth_cnt_sigma2_k <- (sigma2_hat[k]^4 / df_k^4) *
        (mu4_raw - 4 * mu3_raw * mu1_raw +
           6 * mu2_raw * mu1_raw^2 - 3 * mu1_raw^4)

      mu2 <- mu2 + (
        (term1_trace^2 + 6 * term2_trace) * var_sigma2_k -
          (2 * term1_trace * term2_trace + 4 * term3_trace) * third_cnt_sigma2_k +
          (term2_trace^2 + 2 * term4_trace) * fourth_cnt_sigma2_k -
          term2_trace^2 * var_sigma2_k^2
      )
    }

    if (K > 1L) {
      for (k1 in 1:(K - 1L)) {
        for (k2 in (k1 + 1L):K) {
          cap1 <- store_capb[[k1]]
          cap2 <- store_capb[[k2]]

          qnty1 <- cap1 %*% cap2 + cap2 %*% cap1
          qnty2 <- cap1 %*% cap1 %*% cap2 %*% cap2

          qnty1_trace <- sum(diag(qnty1))
          qnty2_trace <- sum(diag(qnty2))

          mu2 <- mu2 + store_var_sigma2[k1] * store_var_sigma2[k2] * (
            qnty1_trace^2 + 4 * qnty2_trace + 2 * sum(diag(qnty1 %*% qnty1))
          )
        }
      }
    }

    myd <- (nrlma - 2 + 2 * nrlma * mu2 / mu1^2) /
      (0.5 * nrlma * mu2 / mu1^2 - 1)
    myc <- (nrlma * myd) / (mu1 * (myd - 2))

    if (!is.finite(myd) || !is.finite(myc) || myd <= 0 || myc <= 0) {
      h0 <- function(parameter) {
        df_para <- exp(parameter[1])
        c_para <- exp(parameter[2])

        (1 / (1 - 2 / df_para) - (c_para * mu1 / nrlma))^2 +
          (
            2 * (1 + nrlma / df_para - 2 / df_para) /
              (nrlma * (1 - 2 / df_para)^2 * (1 - 4 / df_para)) -
              (c_para / nrlma)^2 * mu2
          )^2
      }

      fit_opt <- tryCatch(
        stats::optim(
          c(-0.1, -0.1),
          h0,
          method = "L-BFGS-B",
          lower = c(log(4.01), log(0.001)),
          upper = c(log(100), log(100))
        ),
        error = function(e) NULL
      )

      if (is.null(fit_opt)) {
        pvals[g] <- NA_real_
        next
      }

      myd <- exp(fit_opt$par[1])
      myc <- exp(fit_opt$par[2])

      if (!is.finite(myd) || !is.finite(myc) || myd <= 0 || myc <= 0) {
        pvals[g] <- NA_real_
        next
      }
    }

    p_val <- 1 - stats::pf((myc * stat / nrlma), nrlma, myd)

    if (!is.finite(p_val)) {
      pvals[g] <- NA_real_
      next
    }

    pvals[g] <- max(min(p_val, 1), 0)
  }

  raw_p <- as.numeric(pvals)

  BH_p <- rep(NA_real_, length(raw_p))
  idx <- which(!is.na(raw_p) & is.finite(raw_p) & raw_p >= 0 & raw_p <= 1)

  if (length(idx) > 0) {
    BH_p[idx] <- stats::p.adjust(raw_p[idx], method = "BH")
  }

  qv <- rep(NA_real_, length(raw_p))
  if (length(idx) > 0 && requireNamespace("qvalue", quietly = TRUE)) {
    pv_ok <- raw_p[idx]

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
        qv[idx] <- qobj$qvalues
      }
    }
  }

  out <- data.frame(
    pvalue = raw_p,
    BH = BH_p,
    qvalue = qv,
    row.names = rownames(data_list[[1]])
  )

  return(out)
}
