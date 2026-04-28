#' Test whether rhythm phase differs across multiple conditions.
#'
#' @param data_list List of condition-specific expression matrices.
#' @param timepoint_list List of condition-specific time vectors.
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
#'   The function is applied only to genes that exhibit rhythmicity
#'   in all conditions. A gene is considered rhythmic in a condition if its nominal
#'   TR-test p-value is less than 0.05. Genes that do not meet this criterion
#'   in every condition are assigned \code{NA}.
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
#' # phases for multiple conditions
#' phi1 <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)
#' phi2 <- c(6, 2, 10, 4, 10, 12, 14, 16, 20, 18, 24, 22)
#' phi3 <- c(4, 8, 2, 12, 10, 12, 14, 16, 22, 24, 18, 20)
#'
#' # shared amplitude and mesor
#' amp  <- c(1.2, 1.0, 0.8, 0.6, 0.4, 0.4, 0.2, 0.2, 0.9, 0.7, 0.5, 0.3)
#' mesor <- 1
#'
#' # create simulated expression data
#' data1 <- t(sapply(seq_along(phi1), function(i) {
#'   mesor + amp[i] * cos(2 * pi * (time_vec1 - phi1[i]) / period) +
#'     rnorm(length(time_vec1), 0, 0.2)
#' }))
#'
#' data2 <- t(sapply(seq_along(phi2), function(i) {
#'   mesor + amp[i] * cos(2 * pi * (time_vec2 - phi2[i]) / period) +
#'     rnorm(length(time_vec2), 0, 0.2)
#' }))
#'
#' data3 <- t(sapply(seq_along(phi3), function(i) {
#'   mesor + amp[i] * cos(2 * pi * (time_vec3 - phi3[i]) / period) +
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
#' # combine data and time points as lists
#' data_list <- list(data1, data2, data3)
#' timepoint_list <- list(time_vec1, time_vec2, time_vec3)
#'
#' # run TDP
#' TDP_res <- TDP(data_list = data_list, timepoint_list = timepoint_list, ncond = 3, period = 24)
#'
#' # view results
#' TDP_res
#'
#' # count genes with significant differential phase using nominal p-values below 0.05
#' sum(TDP_res$pvalue < 0.05, na.rm = TRUE)
#'
#' # specify the FDR level
#' fdr_level <- 0.05
#'
#' # count genes with significant differential phase using BH-adjusted p-values
#' sum(TDP_res$BH < fdr_level, na.rm = TRUE)
#'
#' # count genes with significant differential phase using Storey q-values
#' # qvalue is NA if the qvalue package is not installed
#' sum(TDP_res$qvalue < fdr_level, na.rm = TRUE)
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
#' phi1_ns <- c(
#'   2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
#'   22, 24, 2, 4, 6, 8, 10, 12, 14, 16
#' )
#'
#' phi2_ns <- c(
#'   6, 2, 10, 4, 10, 12, 14, 16, 20, 18,
#'   24, 22, 6, 2, 10, 4, 10, 12, 14, 16
#' )
#'
#' phi3_ns <- c(
#'   4, 8, 2, 12, 10, 12, 14, 16, 22, 24,
#'   18, 20, 4, 8, 2, 12, 10, 12, 14, 16
#' )
#'
#' amp_ns <- c(
#'   1.2, 1.0, 0.8, 0.6, 0.4, 0.4, 0.2, 0.2, 0.9, 0.7,
#'   0.5, 0.3, 1.2, 1.0, 0.8, 0.6, 0.4, 0.4, 0.2, 0.2
#' )
#'
#' mesor_ns <- 1
#'
#' data1_ns <- t(sapply(seq_along(phi1_ns), function(i) {
#'   mesor_ns + amp_ns[i] * cos(2 * pi * (time_vec1_ns - phi1_ns[i]) / period) +
#'     rnorm(length(time_vec1_ns), 0, 0.2)
#' }))
#'
#' data2_ns <- t(sapply(seq_along(phi2_ns), function(i) {
#'   mesor_ns + amp_ns[i] * cos(2 * pi * (time_vec2_ns - phi2_ns[i]) / period) +
#'     rnorm(length(time_vec2_ns), 0, 0.2)
#' }))
#'
#' data3_ns <- t(sapply(seq_along(phi3_ns), function(i) {
#'   mesor_ns + amp_ns[i] * cos(2 * pi * (time_vec3_ns - phi3_ns[i]) / period) +
#'     rnorm(length(time_vec3_ns), 0, 0.2)
#' }))
#'
#' rownames(data1_ns) <- paste0("gene", seq_along(phi1_ns))
#' rownames(data2_ns) <- paste0("gene", seq_along(phi2_ns))
#' rownames(data3_ns) <- paste0("gene", seq_along(phi3_ns))
#'
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
#' TDP_res_ns <- TDP(
#'   data_list = data_list_ns,
#'   timepoint_list = timepoint_list_ns,
#'   ncond = 3,
#'   period = 24
#' )
#'
#' TDP_res_ns
#'
#' sum(TDP_res_ns$pvalue < 0.05, na.rm = TRUE)
#'
#' fdr_level <- 0.05
#'
#' sum(TDP_res_ns$BH < fdr_level, na.rm = TRUE)
#'
#' sum(TDP_res_ns$qvalue < fdr_level, na.rm = TRUE)
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
#' phi1_miss <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)
#' phi2_miss <- c(6, 2, 10, 4, 10, 12, 14, 16, 20, 18, 24, 22)
#' phi3_miss <- c(4, 8, 2, 12, 10, 12, 14, 16, 22, 24, 18, 20)
#'
#' amp_miss <- c(1.2, 1.0, 0.8, 0.6, 0.4, 0.4, 0.2, 0.2, 0.9, 0.7, 0.5, 0.3)
#' mesor_miss <- 1
#'
#' data1_miss <- t(sapply(seq_along(phi1_miss), function(i) {
#'   mesor_miss + amp_miss[i] * cos(2 * pi * (time_vec1_miss - phi1_miss[i]) / period) +
#'     rnorm(length(time_vec1_miss), 0, 0.2)
#' }))
#'
#' data2_miss <- t(sapply(seq_along(phi2_miss), function(i) {
#'   mesor_miss + amp_miss[i] * cos(2 * pi * (time_vec2_miss - phi2_miss[i]) / period) +
#'     rnorm(length(time_vec2_miss), 0, 0.2)
#' }))
#'
#' data3_miss <- t(sapply(seq_along(phi3_miss), function(i) {
#'   mesor_miss + amp_miss[i] * cos(2 * pi * (time_vec3_miss - phi3_miss[i]) / period) +
#'     rnorm(length(time_vec3_miss), 0, 0.2)
#' }))
#'
#' rownames(data1_miss) <- paste0("gene", seq_along(amp_miss))
#' rownames(data2_miss) <- paste0("gene", seq_along(amp_miss))
#' rownames(data3_miss) <- paste0("gene", seq_along(amp_miss))
#'
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
#' # For gene11 and gene12 in condition 1, keep only one observation
#' # at 2 h and 6 h; set the remaining observations at those time points to NA.
#' # time 2 h -> columns 1, 2, 3
#' # time 6 h -> columns 4, 5, 6
#' data1_miss[c("gene11", "gene12"), c(2, 3, 5, 6)] <- NA
#'
#' data_list_miss <- list(data1_miss, data2_miss, data3_miss)
#' timepoint_list_miss <- list(time_vec1_miss, time_vec2_miss, time_vec3_miss)
#'
#' TDP_res_miss <- TDP(data_list = data_list_miss,
#'                     timepoint_list = timepoint_list_miss,
#'                     ncond = 3, period = 24)
#' TDP_res_miss
#'
#' sum(TDP_res_miss$pvalue < 0.05, na.rm = TRUE)
#'
#' fdr_level <- 0.05
#'
#' sum(TDP_res_miss$BH < fdr_level, na.rm = TRUE)
#'
#' sum(TDP_res_miss$qvalue < fdr_level, na.rm = TRUE)
#'
#' @export

TDP <- function(data_list, timepoint_list, ncond, period = 24){

  ncond_data <- length(data_list)

  if (!missing(ncond) && ncond != ncond_data) {
    warning("'ncond' does not match length(data_list); using length(data_list).")
  }
  ncond <- ncond_data

  if (length(timepoint_list) != ncond) {
    stop(paste("Expected", ncond, "sets of timepoints, but got", length(timepoint_list)))
  }

  data_list <- lapply(data_list, as.matrix)

  check_result <- mapply(function(data, time) {
    ncol(data) == length(time)
  }, data_list, timepoint_list)

  if (any(!check_result)) {
    stop("Column count mismatch: at least one dataset does not match its timepoint vector.")
  }

  row_counts <- sapply(data_list, nrow)
  if (length(unique(row_counts)) != 1) {
    stop("All datasets in `data_list` must have the same number of rows.")
  }

  p <- unique(row_counts)
  results <- matrix(NA_real_, nrow = p, ncol = 4)

  alpha_tr <- 0.05
  tr_pvalues_list <- lapply(seq_len(ncond), function(i) {
    tr_res <- TR(
      data_list = data_list[[i]],
      time_vec  = timepoint_list[[i]],
      period    = period
    )

    if (is.data.frame(tr_res)) {
      if (!"pvalue" %in% colnames(tr_res)) {
        stop("TR() returned a data.frame but no `pvalue` column was found.")
      }
      as.numeric(tr_res$pvalue)
    } else if (is.numeric(tr_res)) {
      as.numeric(tr_res)
    } else {
      stop("TR() must return either a numeric vector or a data.frame with column `pvalue`.")
    }
  })

  rhythmic_all <- Reduce(
    `&`,
    lapply(tr_pvalues_list, function(pv) !is.na(pv) & is.finite(pv) & pv < alpha_tr)
  )

  Lmat <- matrix(0, ncol = ncond, nrow = (ncond - 1), byrow = TRUE)
  Lmat[, 1] <- 1
  for (r in 1:(ncond - 1)) {
    Lmat[r, r + 1] <- -1
  }
  nrlma <- nrow(Lmat)

  row_nm <- rownames(data_list[[1]])
  if (is.null(row_nm)) {
    row_nm <- paste0("gene", seq_len(p))
  }

  for (i in seq_len(p)) {

    if (!rhythmic_all[i]) {
      next
    }

    theta_vec <- rep(NA_real_, ncond)
    var_error_vec <- rep(NA_real_, ncond)
    deg_frdm_vec <- rep(NA_real_, ncond)
    amat <- matrix(0, ncol = ncond, nrow = ncond)
    varcov <- matrix(0, ncol = ncond, nrow = ncond)
    xtxinv_list_gene <- vector("list", ncond)
    capd_list <- vector("list", ncond)

    gene_ok <- TRUE

    for (myk in seq_len(ncond)) {

      y_all <- as.numeric(data_list[[myk]][i, ])
      t_all <- as.numeric(timepoint_list[[myk]])

      ok <- is.finite(y_all) & is.finite(t_all)
      y <- y_all[ok]
      t <- t_all[ok]

      n_obs <- length(y)
      if (n_obs <= 3) {
        gene_ok <- FALSE
        break
      }

      mymat <- cbind(
        1,
        cos(2 * pi * t / period),
        sin(2 * pi * t / period)
      )

      if (qr(mymat)$rank < 3) {
        gene_ok <- FALSE
        break
      }

      xtx <- crossprod(mymat)
      xtxinv <- tryCatch(solve(xtx), error = function(e) NULL)
      if (is.null(xtxinv)) {
        gene_ok <- FALSE
        break
      }

      gamma_hat <- xtxinv %*% crossprod(mymat, y)
      resid <- as.numeric(y - mymat %*% gamma_hat)

      deg_frdm <- nrow(mymat) - ncol(mymat)
      if (deg_frdm <= 0) {
        gene_ok <- FALSE
        break
      }

      var_error <- sum(resid^2) / deg_frdm
      if (!is.finite(var_error) || var_error < 0) {
        gene_ok <- FALSE
        break
      }

      denom <- as.numeric(gamma_hat[2]^2 + gamma_hat[3]^2)
      if (!is.finite(denom) || denom <= .Machine$double.eps) {
        gene_ok <- FALSE
        break
      }

      theta_vec[myk] <- atan2(gamma_hat[3], gamma_hat[2])
      capd_list[[myk]] <- c(
        0,
        -as.numeric(gamma_hat[3]) / denom,
        as.numeric(gamma_hat[2]) / denom
      )

      amat[myk, myk] <- as.numeric(
        t(capd_list[[myk]]) %*% xtxinv %*% capd_list[[myk]]
      )
      varcov[myk, myk] <- var_error * amat[myk, myk]

      xtxinv_list_gene[[myk]] <- xtxinv
      deg_frdm_vec[myk] <- deg_frdm
      var_error_vec[myk] <- var_error
    }

    if (!gene_ok) {
      next
    }

    b4inv <- Lmat %*% varcov %*% t(Lmat)
    if (any(!is.finite(b4inv))) {
      next
    }

    Lcov <- tryCatch(solve(b4inv), error = function(e) NULL)
    if (is.null(Lcov)) {
      next
    }

    diff_theta <- Lmat %*% theta_vec
    newtest100 <- as.numeric(t(diff_theta) %*% Lcov %*% diff_theta)
    if (!is.finite(newtest100) || newtest100 < 0) {
      next
    }

    mu1 <- nrlma
    mu2 <- 2 * nrlma
    dummy_varcov <- matrix(0, ncol = ncond, nrow = ncond)
    store_var_sigma2_kg <- rep(NA_real_, ncond)
    store_capg_list <- vector("list", ncond)

    for (k1 in seq_len(ncond)) {

      deg_frdm <- deg_frdm_vec[k1]
      dummy_varcov[k1, k1] <- amat[k1, k1]

      capgk1 <- Lmat %*% dummy_varcov %*% t(Lmat) %*% Lcov
      term2 <- capgk1 %*% capgk1
      term3 <- term2 %*% capgk1
      term4 <- term3 %*% capgk1

      term1_trace <- sum(diag(capgk1))
      term2_trace <- sum(diag(term2))
      term3_trace <- sum(diag(term3))
      term4_trace <- sum(diag(term4))

      store_capg_list[[k1]] <- capgk1

      mu1 <- mu1 + term2_trace * 2 * var_error_vec[k1]^2 / deg_frdm

      dummy_varcov[k1, k1] <- 0

      mu1_raw <- deg_frdm
      mu2_raw <- 2 * (deg_frdm / 2 + 1) * deg_frdm
      mu3_raw <- 4 * (deg_frdm / 2 + 2) * (deg_frdm / 2 + 1) * deg_frdm
      mu4_raw <- 8 * (deg_frdm / 2 + 3) * (deg_frdm / 2 + 2) * (deg_frdm / 2 + 1) * deg_frdm

      var_sigma2_kg <- 2 * var_error_vec[k1]^2 / deg_frdm
      store_var_sigma2_kg[k1] <- var_sigma2_kg

      thrd_cnt_mnt_sigma2kg <- ((var_error_vec[k1]^3 / deg_frdm^3) *
                                  (mu3_raw - 3 * mu2_raw * mu1_raw + 2 * mu1_raw^3))

      frth_cnt_mnt_sigma2kg <- ((var_error_vec[k1]^4 / deg_frdm^4) *
                                  (mu4_raw - 4 * mu3_raw * mu1_raw + 6 * mu2_raw * mu1_raw^2 - 3 * mu1_raw^4))

      mu2 <- mu2 + (
        (term1_trace * term1_trace + 6 * term2_trace) * var_sigma2_kg
        - (2 * term1_trace * term2_trace + 4 * term3_trace) * thrd_cnt_mnt_sigma2kg
        + (term2_trace * term2_trace + 2 * term4_trace) * frth_cnt_mnt_sigma2kg
        - term2_trace * term2_trace * var_sigma2_kg^2
      )
    }

    if (!is.finite(mu1) || !is.finite(mu2) || mu1 == 0) {
      next
    }

    for (k1 in 1:(ncond - 1)) {
      for (k1new in (k1 + 1):ncond) {
        qnty1 <- store_capg_list[[k1]] %*% store_capg_list[[k1new]] +
          store_capg_list[[k1new]] %*% store_capg_list[[k1]]
        qnty2 <- store_capg_list[[k1]] %*% store_capg_list[[k1]] %*%
          store_capg_list[[k1new]] %*% store_capg_list[[k1new]]

        qnty1_trace <- sum(diag(qnty1))
        qnty2_trace <- sum(diag(qnty2))

        mu2 <- mu2 + store_var_sigma2_kg[k1] * store_var_sigma2_kg[k1new] * (
          qnty1_trace * qnty1_trace +
            4 * qnty2_trace +
            2 * sum(diag(qnty1 %*% qnty1))
        )
      }
    }

    if (!is.finite(mu2)) {
      next
    }

    myd <- (nrlma - 2 + 2 * nrlma * mu2 / mu1^2) /
      (0.5 * nrlma * mu2 / mu1^2 - 1)
    myc <- (nrlma * myd) / (mu1 * (myd - 2))

    if (!is.finite(myd) || !is.finite(myc) || myd < 0 || myc < 0) {

      h0 <- function(parameter) {
        df.para <- exp(parameter[1])
        c.para <- exp(parameter[2])

        (1 / (1 - 2 / df.para) - (c.para * mu1 / nrlma))^2 +
          (
            2 * (1 + nrlma / df.para - 2 / df.para) /
              (nrlma * (1 - 2 / df.para)^2 * (1 - 4 / df.para)) -
              (c.para / nrlma)^2 * mu2
          )^2
      }

      outh0 <- tryCatch(
        optim(
          c(-0.1, -0.1), h0, method = "L-BFGS-B",
          lower = c(log(4.01), log(0.001)),
          upper = c(log(100), log(100))
        ),
        error = function(e) NULL
      )

      if (is.null(outh0)) {
        next
      }

      myd <- exp(outh0$par[1])
      myc <- exp(outh0$par[2])
    }

    if (!is.finite(myd) || !is.finite(myc) || myd <= 0 || myc <= 0) {
      next
    }

    p_val <- 1 - pf((myc * newtest100 / nrlma), nrlma, myd)
    if (!is.finite(p_val)) {
      next
    }

    results[i, ] <- c(newtest100, nrlma, myd, p_val)
  }

  raw_p <- results[, 4]

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
    row.names = row_nm
  )

  return(out)
}


