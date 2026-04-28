#' Test whether rhythm amplitude differs across multiple conditions.
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
#' # amplitudes for multiple conditions
#' amp1 <- c(1.2, 1.0, 0.8, 0.6, 0.4, 0.4, 0.2, 0.2, 1.1, 0.7, 0.5, 0.3)
#' amp2 <- c(0.6, 1.4, 0.3, 1.0, 0.4, 0.4, 0.2, 0.2, 0.9, 1.2, 0.6, 0.4)
#' amp3 <- c(1.5, 0.5, 1.1, 0.2, 0.4, 0.4, 0.2, 0.2, 1.3, 0.8, 0.7, 0.3)
#'
#' # shared phase and mesor
#' phi  <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)
#' mesor <- 1
#'
#' # create simulated expression data
#' data1 <- t(sapply(seq_along(amp1), function(i) {
#'   mesor + amp1[i] * cos(2 * pi * (time_vec1 - phi[i]) / period) +
#'     rnorm(length(time_vec1), 0, 0.2)
#' }))
#'
#' data2 <- t(sapply(seq_along(amp2), function(i) {
#'   mesor + amp2[i] * cos(2 * pi * (time_vec2 - phi[i]) / period) +
#'     rnorm(length(time_vec2), 0, 0.2)
#' }))
#'
#' data3 <- t(sapply(seq_along(amp3), function(i) {
#'   mesor + amp3[i] * cos(2 * pi * (time_vec3 - phi[i]) / period) +
#'     rnorm(length(time_vec3), 0, 0.2)
#' }))
#'
#' rownames(data1) <- paste0("gene", seq_along(amp1))
#' rownames(data2) <- paste0("gene", seq_along(amp2))
#' rownames(data3) <- paste0("gene", seq_along(amp3))
#' colnames(data1) <- paste0("Time ", time_vec1, ", rep ", rep(1:3, times = 6))
#' colnames(data2) <- paste0("Time ", time_vec2, ", rep ", rep(1:3, times = 6))
#' colnames(data3) <- paste0("Time ", time_vec3, ", rep ", rep(1:3, times = 6))
#'
#' # combine data and time points as lists
#' data_list <- list(data1, data2, data3)
#' timepoint_list <- list(time_vec1, time_vec2, time_vec3)
#'
#' # run TDA
#' TDA_res <- TDA(data_list = data_list, timepoint_list = timepoint_list, ncond = 3, period = 24)
#'
#' # view results
#' TDA_res
#'
#' # count genes with significant differential amplitude using nominal p-values below 0.05
#' sum(TDA_res$pvalue < 0.05, na.rm = TRUE)
#'
#' # specify the FDR level
#' fdr_level <- 0.05
#'
#' # count genes with significant differential amplitude using BH-adjusted p-values
#' sum(TDA_res$BH < fdr_level, na.rm = TRUE)
#'
#' # count genes with significant differential amplitude using Storey q-values
#' # qvalue is NA if the qvalue package is not installed
#' sum(TDA_res$qvalue < fdr_level, na.rm = TRUE)
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
#' amp1_ns <- c(
#'   1.2, 1.0, 0.8, 0.6, 0.4, 0.4, 0.2, 0.2, 1.1, 0.7,
#'   0.5, 0.3, 1.2, 1.0, 0.8, 0.6, 0.4, 0.4, 0.2, 0.2
#' )
#'
#' amp2_ns <- c(
#'   0.6, 1.4, 0.3, 1.0, 0.4, 0.4, 0.2, 0.2, 0.9, 1.2,
#'   0.6, 0.4, 0.6, 1.4, 0.3, 1.0, 0.4, 0.4, 0.2, 0.2
#' )
#'
#' amp3_ns <- c(
#'   1.5, 0.5, 1.1, 0.2, 0.4, 0.4, 0.2, 0.2, 1.3, 0.8,
#'   0.7, 0.3, 1.5, 0.5, 1.1, 0.2, 0.4, 0.4, 0.2, 0.2
#' )
#'
#' phi_ns <- c(
#'   2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
#'   22, 24, 2, 4, 6, 8, 10, 12, 14, 16
#' )
#'
#' mesor_ns <- 1
#'
#' data1_ns <- t(sapply(seq_along(amp1_ns), function(i) {
#'   mesor_ns + amp1_ns[i] * cos(2 * pi * (time_vec1_ns - phi_ns[i]) / period) +
#'     rnorm(length(time_vec1_ns), 0, 0.2)
#' }))
#'
#' data2_ns <- t(sapply(seq_along(amp2_ns), function(i) {
#'   mesor_ns + amp2_ns[i] * cos(2 * pi * (time_vec2_ns - phi_ns[i]) / period) +
#'     rnorm(length(time_vec2_ns), 0, 0.2)
#' }))
#'
#' data3_ns <- t(sapply(seq_along(amp3_ns), function(i) {
#'   mesor_ns + amp3_ns[i] * cos(2 * pi * (time_vec3_ns - phi_ns[i]) / period) +
#'     rnorm(length(time_vec3_ns), 0, 0.2)
#' }))
#'
#' rownames(data1_ns) <- paste0("gene", seq_along(amp1_ns))
#' rownames(data2_ns) <- paste0("gene", seq_along(amp1_ns))
#' rownames(data3_ns) <- paste0("gene", seq_along(amp1_ns))
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
#' TDA_res_ns <- TDA(
#'   data_list = data_list_ns,
#'   timepoint_list = timepoint_list_ns,
#'   ncond = 3,
#'   period = 24
#' )
#'
#' TDA_res_ns
#'
#' sum(TDA_res_ns$pvalue < 0.05, na.rm = TRUE)
#'
#' fdr_level <- 0.05
#'
#' sum(TDA_res_ns$BH < fdr_level, na.rm = TRUE)
#'
#' sum(TDA_res_ns$qvalue < fdr_level, na.rm = TRUE)
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
#' amp1_miss <- c(
#'   1.2, 1.0, 0.8, 0.6, 0.4, 0.4, 0.2, 0.2, 1.1, 0.7,
#'   0.5, 0.3, 1.2, 1.0, 0.8, 0.6, 0.4, 0.4, 0.2, 0.2
#' )
#'
#' amp2_miss <- c(
#'   0.6, 1.4, 0.3, 1.0, 0.4, 0.4, 0.2, 0.2, 0.9, 1.2,
#'   0.6, 0.4, 0.6, 1.4, 0.3, 1.0, 0.4, 0.4, 0.2, 0.2
#' )
#'
#' amp3_miss <- c(
#'   1.5, 0.5, 1.1, 0.2, 0.4, 0.4, 0.2, 0.2, 1.3, 0.8,
#'   0.7, 0.3, 1.5, 0.5, 1.1, 0.2, 0.4, 0.4, 0.2, 0.2
#' )
#'
#' phi_miss <- c(
#'   2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
#'   22, 24, 2, 4, 6, 8, 10, 12, 14, 16
#' )
#'
#' mesor_miss <- 1
#'
#' data1_miss <- t(sapply(seq_along(amp1_miss), function(i) {
#'   mesor_miss + amp1_miss[i] * cos(2 * pi * (time_vec1_miss - phi_miss[i]) / period) +
#'     rnorm(length(time_vec1_miss), 0, 0.2)
#' }))
#'
#' data2_miss <- t(sapply(seq_along(amp2_miss), function(i) {
#'   mesor_miss + amp2_miss[i] * cos(2 * pi * (time_vec2_miss - phi_miss[i]) / period) +
#'     rnorm(length(time_vec2_miss), 0, 0.2)
#' }))
#'
#' data3_miss <- t(sapply(seq_along(amp3_miss), function(i) {
#'   mesor_miss + amp3_miss[i] * cos(2 * pi * (time_vec3_miss - phi_miss[i]) / period) +
#'     rnorm(length(time_vec3_miss), 0, 0.2)
#' }))
#'
#' rownames(data1_miss) <- paste0("gene", seq_along(amp1_miss))
#' rownames(data2_miss) <- paste0("gene", seq_along(amp1_miss))
#' rownames(data3_miss) <- paste0("gene", seq_along(amp1_miss))
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
#' # For gene19 and gene20 in condition 1, keep only one observation
#' # at 2 h and 6 h; set the remaining observations at those time points to NA.
#' # time 2 h -> columns 1, 2, 3
#' # time 6 h -> columns 4, 5, 6
#' data1_miss[c("gene19", "gene20"), c(2, 3, 5, 6)] <- NA
#'
#' data_list_miss <- list(data1_miss, data2_miss, data3_miss)
#' timepoint_list_miss <- list(time_vec1_miss, time_vec2_miss, time_vec3_miss)
#'
#' TDA_res_miss <- TDA(
#'   data_list = data_list_miss,
#'   timepoint_list = timepoint_list_miss,
#'   ncond = 3,
#'   period = 24
#' )
#'
#' TDA_res_miss
#'
#' sum(TDA_res_miss$pvalue < 0.05, na.rm = TRUE)
#'
#' fdr_level <- 0.05
#'
#' sum(TDA_res_miss$BH < fdr_level, na.rm = TRUE)
#'
#' sum(TDA_res_miss$qvalue < fdr_level, na.rm = TRUE)
#'
#' @export

TDA <- function(data_list, timepoint_list, ncond, period = 24){

  ncond <- length(data_list)

  if (length(data_list) != ncond) {
    stop(paste("Expected", ncond, "datasets, but got", length(data_list)))
  }
  if (length(timepoint_list) != ncond) {
    stop(paste("Expected", ncond, "set of timepoints, but got", length(timepoint_list)))
  }


  check_result <- mapply(function(data, time) {
    ncol(data) == length(time)
  }, data_list, timepoint_list)



  if (any(!check_result)) {
    stop("Column count mismatch: at least one dataset does not match the expected column count.")
  }

  row_counts <- sapply(data_list, nrow)
  if (length(unique(row_counts)) != 1) {
    stop("All datasets in `data_list` must have the same number of rows.")
  }

  p <- unique(row_counts)
  results <- matrix(NA_real_, nrow = p, ncol = 4)

  alpha_tr = 0.05
  tr_pvalues_list <- NULL

  if (is.null(tr_pvalues_list)) {
    tr_pvalues_list <- lapply(seq_len(ncond), function(i) {
      tr_res <- TR(
        data_list = data_list[[i]],
        time_vec  = timepoint_list[[i]],
        period = period
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
  }

  rhythmic_all <- Reduce(
    `&`,
    lapply(tr_pvalues_list, function(pv) !is.na(pv) & pv < alpha_tr)
  )

  Lmat <- matrix(0, ncol = ncond, nrow = (ncond - 1), byrow = TRUE)
  Lmat[, 1] <- 1
  for (r in 1:(ncond - 1)) {
    Lmat[r, (r + 1)] <- -1
  }
  nrlma <- nrow(Lmat)

  mymat_list <- vector("list", ncond)
  mymat_nr_list <- vector("list", ncond)
  mymat_nc_list <- vector("list", ncond)
  xtxinv_list <- vector("list", ncond)
  xtxinvxt_list <- vector("list", ncond)

  for (i in seq_len(p)) {

    if (!rhythmic_all[i]) {
      results[i, 4] <- NA_real_
      next
    }


  for(myk in 1:ncond){
    nobseachcond=length(timepoint_list[[myk]])
    mymat=matrix(0, nrow=nobseachcond, ncol= 3)

    mymat[, 1]=1

    mymat[, 2]=cos(2*pi*timepoint_list[[myk]]/24)
    mymat[, 3]=sin(2*pi*timepoint_list[[myk]]/24)
    ncmymat=ncol(mymat)
    nrmymat=nrow(mymat)
    xtx=t(mymat) %*% mymat
    xtxinv=solve(xtx)
    xtxinvxt=xtxinv %*% t(mymat)

    mymat_list[[myk]] <- mymat
    mymat_nr_list[[myk]] <- nrmymat
    mymat_nc_list[[myk]] <- ncmymat
    xtxinv_list[[myk]] <- xtxinv
    xtxinvxt_list[[myk]] <- xtxinvxt
  }

  gamma_vec=rep(0, (ncond*3));
  var_error_vec=rep(0, ncond);
  amat=varcov=matrix(0, ncol=ncond, nrow=ncond)


    theta_vec<- rep(0, ncond)
    capd_list<- vector("list", ncond)

    for(myk in 1:ncond){
      newy <- data_list[[myk]][i, ]
      gamma_hat <- xtxinvxt_list[[myk]] %*% as.vector(t(newy));
      var_error=(t(as.vector(t(newy))-mymat_list[[myk]]%*%gamma_hat)%*%(
        as.vector(t(newy))-mymat_list[[myk]]%*%gamma_hat)/(mymat_nr_list[[myk]]-mymat_nc_list[[myk]]))
      k1= (myk-1)*3+1;
      k2=myk*3;

      theta_vec[myk]=(sum(gamma_hat[-1]^2))^(1/3)
      capd_list[[myk]]<- (2/3)*c(0, gamma_hat[2]/theta_vec[myk]^2, gamma_hat[3]/theta_vec[myk]^2)

      amat[myk, myk]= (t(capd_list[[myk]])%*%
                         xtxinv_list[[myk]]%*%capd_list[[myk]])

      varcov[myk, myk]= as.numeric(var_error)* amat[myk, myk]

      gamma_vec[k1:k2]=gamma_hat
      var_error_vec[myk]=var_error

    }

    b4inv=Lmat %*% varcov %*% t(Lmat)

    Lcov <- solve(b4inv);

    newtest100 <-(as.numeric(t(Lmat %*% theta_vec)%*%Lcov%*%(Lmat %*% theta_vec)))

    mu1= nrlma
    mu2= 2*nrlma
    dummy_varcov=0*varcov
    store_var_sigma2_kg= rep(0, ncond)
    store_capg_list <- vector("list", ncond)

    for(k1 in 1:ncond){

      deg_frdm= (mymat_nr_list[[k1]]-mymat_nc_list[[k1]])
      k2=3*(k1-1)+1
      k3=3*k1
      dummy_varcov[k1, k1]=amat[k1, k1]
      capgk1= Lmat %*% dummy_varcov%*%t(Lmat)%*%Lcov
      term2= capgk1%*%capgk1
      term3= term2%*%capgk1
      term4= term3%*%capgk1
      term1_trace= sum(diag(capgk1))
      term2_trace= sum(diag(term2))
      term3_trace= sum(diag(term3))
      term4_trace= sum(diag(term4))

      store_capg_list[[k1]] <- capgk1

      mu1=mu1+ (term2_trace)*2*var_error_vec[k1]^2/deg_frdm

      dummy_varcov[k1, k1]=0

      mu1_raw= deg_frdm
      mu2_raw= 2*(deg_frdm/2+1)*deg_frdm
      mu3_raw= 4*(deg_frdm/2+2)*(deg_frdm/2+1)*deg_frdm
      mu4_raw= 8*(deg_frdm/2+3)*(deg_frdm/2+2)*(deg_frdm/2+1)*deg_frdm
      var_sigma2_kg= 2*var_error_vec[k1]^2/deg_frdm
      store_var_sigma2_kg[k1]=var_sigma2_kg


      thrd_cnt_mnt_sigma2kg=((var_error_vec[k1]^3/deg_frdm^3)*
                               (mu3_raw-3*mu2_raw*mu1_raw+2*mu1_raw^3))
      frth_cnt_mnt_sigma2kg=((var_error_vec[k1]^4/deg_frdm^4)*
                               (mu4_raw-4*mu3_raw*mu1_raw+6*mu2_raw*mu1_raw^2-3*mu1_raw^4))


      mu2=mu2+ (
        (term1_trace*term1_trace+6*term2_trace)*var_sigma2_kg
        -(2*term1_trace*term2_trace+4*term3_trace)*thrd_cnt_mnt_sigma2kg
        +(term2_trace*term2_trace+2*term4_trace)*frth_cnt_mnt_sigma2kg
        - term2_trace*term2_trace*var_sigma2_kg^2
      )

    }

    for(k1 in 1:(ncond-1)){
      for(k1new in (k1+1):ncond){
        qnty1=( store_capg_list[[k1]]%*%store_capg_list[[k1new]]+
                  store_capg_list[[k1new]]%*%store_capg_list[[k1]])
        qnty2=(store_capg_list[[k1]]%*%store_capg_list[[k1]]%*%
                 store_capg_list[[k1new]]%*%store_capg_list[[k1new]])
        qnty1_trace= sum(diag(qnty1))
        qnty2_trace= sum(diag(qnty2))
        mu2= mu2+store_var_sigma2_kg[k1]*store_var_sigma2_kg[k1new]*(
          qnty1_trace*qnty1_trace+4*qnty2_trace+ 2*sum(diag(qnty1%*%qnty1)))
      }
    }

    myd= (nrlma-2+2*nrlma*mu2/mu1^2)/(0.5*nrlma*mu2/mu1^2-1)
    myc= (nrlma*myd)/(mu1*(myd-2))
    if(myd<0 |myc<0){

      h0=function(parameter){
        df.para=exp(parameter[1])
        c.para=exp(parameter[2])
        return((1/(1-2/df.para)- (c.para*mu1/nrlma))^2+
                 (2*(1+nrlma/df.para-2/df.para)/(nrlma*(1-2/df.para)^2*(1-4/df.para))
                  -(c.para/nrlma)^2*mu2)^2
        )

      }
      outh0=optim(c(-0.1, -0.1), h0, method="L-BFGS-B",

                  lower=c(log(4.01), log(0.001)), upper=c(log(100), log(100)))

      myd=exp(outh0$par[1])
      myc=exp(outh0$par[2])
    }

    p_val=1-pf((myc*newtest100/nrlma), nrlma, myd)
    results[i, ]=c(newtest100,  nrlma, myd, p_val)
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
    row.names = rownames(data_list[[1]])
  )

  return(out)
}
