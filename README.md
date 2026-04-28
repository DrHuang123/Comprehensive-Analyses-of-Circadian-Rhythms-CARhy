# CARhy

<!-- badges: start -->
<!-- badges: end -->

**CARhy** (*Comprehensive Analyses of Circadian Rhythms*) is an R package for analyzing circadian rhythms in transcriptomic experiments with multiple experimental conditions. It provides functions for rhythmicity testing, differential rhythmicity testing, differential mesor/amplitude/phase testing, rhythm-parameter estimation, and multiple-testing correction.

## Features

- Test rhythmicity within a single condition.
- Test differential rhythmicity across multiple conditions.
- Test differential rhythm amplitude, phase, and mesor across multiple conditions.
- Estimate mesor, amplitude, phase, standard errors, and confidence intervals.
- Supports balanced designs, unevenly spaced sampling times, unequal numbers of replicates, and missing values.
- Reports raw p-values, Benjamini-Hochberg adjusted p-values, and Storey q-values when the **qvalue** package is available.
- Optional RNA-seq count preprocessing using low-expression filtering and TMM normalization via **edgeR**.

## Installation

Install the CARhy (Comprehensive Analyses of Circadian Rhythms):

```r
install.packages("remotes")
remotes::install_github("DrHuang123/Comprehensive-Analyses-of-Circadian-Rhythms-CARhy")
```

For optional RNA-seq preprocessing and Storey q-value estimation, install the Bioconductor dependencies:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("edgeR", "qvalue"))
```

## Quick start

```r
library(CARhy)

set.seed(1)
period <- 24

# Six sampling times with three replicates at each time point
time_vec <- rep(c(2, 6, 10, 14, 18, 22), each = 3)

# Simulate expression data for eight genes
amp <- c(1.2, 1.0, 0.8, 0.6, 0, 0, 0.4, 0.2)
phi <- c(2, 4, 6, 8, 0, 0, 10, 12)
mesor <- 1

expr <- t(sapply(seq_along(amp), function(i) {
  mesor + amp[i] * cos(2 * pi * (time_vec - phi[i]) / period) +
    rnorm(length(time_vec), 0, 0.2)
}))
rownames(expr) <- paste0("gene", seq_along(amp))
colnames(expr) <- paste0("Time ", time_vec, ", rep ", rep(1:3, times = 6))

# Test rhythmicity within one condition
TR_res <- TR(expr, time_vec, period = 24)
head(TR_res)

# Count significant rhythmic genes using BH-adjusted p-values
sum(TR_res$BH < 0.05, na.rm = TRUE)

# Estimate rhythm parameters
param_res <- params_output(expr, time_vec, period = 24)
head(param_res)
```

## Main functions

| Function | Purpose | Main input | Output |
| --- | --- | --- | --- |
| `preprocess_counts()` | Preprocess RNA-seq count data before rhythm analysis | Raw count matrix | Normalized expression matrix or `edgeR::DGEList` |
| `TR()` | Test rhythmicity within one condition | Gene-by-sample expression matrix and sampling times | `pvalue`, `BH`, `qvalue` |
| `TDR()` | Test whether rhythmicity differs across multiple conditions | List of expression matrices and time vectors | `pvalue`, `BH`, `qvalue` |
| `TDA()` | Test whether rhythm amplitude differs across multiple conditions | List of expression matrices and time vectors | `pvalue`, `BH`, `qvalue` |
| `TDP()` | Test whether rhythm phase differs across multiple conditions | List of expression matrices and time vectors | `pvalue`, `BH`, `qvalue` |
| `TDM()` | Test whether mesor differs across multiple conditions | List of expression matrices and time vectors | `pvalue`, `BH`, `qvalue` |
| `params_output()` | Fit the rhythm model and extract parameters | Expression matrix and sampling times | Mesor, amplitude, phase, SEs, and confidence intervals |

## RNA-seq count preprocessing

`preprocess_counts()` provides an optional workflow for RNA-seq count data: low-expression filtering with `edgeR::filterByExpr()`, followed by TMM normalization with **edgeR**.

```r
set.seed(123)
counts <- matrix(rnbinom(2000, mu = 50, size = 1), nrow = 200, ncol = 10)
group <- rep(c("A", "B"), each = 5)

expr <- preprocess_counts(
  counts,
  filter = TRUE,
  group = group,
  method = "TMM",
  output = "logCPM"
)

head(expr)
```

You may also use any other suitable normalization workflow before applying CARhy functions.

## Multi-condition example

```r
set.seed(1)
period <- 24

time_vec1 <- rep(c(2, 6, 10, 14, 18, 22), each = 3)
time_vec2 <- rep(c(2, 6, 10, 14, 18, 22), each = 3)
time_vec3 <- rep(c(2, 6, 10, 14, 18, 22), each = 3)

amp1 <- c(1.2, 1.0, 0.8, 0.6, 0.4, 0.4, 0.2, 0.2)
amp2 <- c(0.6, 1.4, 0.3, 1.0, 0.4, 0.4, 0.2, 0.2)
amp3 <- c(1.5, 0.5, 1.1, 0.2, 0.4, 0.4, 0.2, 0.2)

phi1 <- c(2, 4, 6, 8, 10, 12, 14, 16)
phi2 <- c(6, 2, 10, 4, 10, 12, 14, 16)
phi3 <- c(4, 8, 2, 12, 10, 12, 14, 16)

mesor1 <- rep(1.0, 8)
mesor2 <- rep(1.3, 8)
mesor3 <- rep(0.8, 8)

data1 <- t(sapply(seq_along(amp1), function(i) {
  mesor1[i] + amp1[i] * cos(2 * pi * (time_vec1 - phi1[i]) / period) +
    rnorm(length(time_vec1), 0, 0.2)
}))
data2 <- t(sapply(seq_along(amp2), function(i) {
  mesor2[i] + amp2[i] * cos(2 * pi * (time_vec2 - phi2[i]) / period) +
    rnorm(length(time_vec2), 0, 0.2)
}))
data3 <- t(sapply(seq_along(amp3), function(i) {
  mesor3[i] + amp3[i] * cos(2 * pi * (time_vec3 - phi3[i]) / period) +
    rnorm(length(time_vec3), 0, 0.2)
}))

rownames(data1) <- rownames(data2) <- rownames(data3) <- paste0("gene", seq_along(amp1))

data_list <- list(data1, data2, data3)
timepoint_list <- list(time_vec1, time_vec2, time_vec3)

# Differential rhythmicity
TDR_res <- TDR(data_list, timepoint_list, ncond = 3, period = 24)

# Differential amplitude, phase, and mesor
TDA_res <- TDA(data_list, timepoint_list, ncond = 3, period = 24)
TDP_res <- TDP(data_list, timepoint_list, ncond = 3, period = 24)
TDM_res <- TDM(data_list, timepoint_list, ncond = 3, period = 24)

sum(TDR_res$BH < 0.05, na.rm = TRUE)
sum(TDA_res$BH < 0.05, na.rm = TRUE)
sum(TDP_res$BH < 0.05, na.rm = TRUE)
sum(TDM_res$BH < 0.05, na.rm = TRUE)
```

## Rhythm-parameter estimation

`params_output()` fits the model

```text
expression ~ cos(2 * pi * time / period) + sin(2 * pi * time / period)
```

for each gene and returns mesor, amplitude, phase, standard errors, and confidence intervals.

```r
param_res <- params_output(
  expr_mat = expr,
  time_vec = time_vec,
  period = 24,
  conf_level = 0.95
)

head(param_res)
```

Missing values are removed together with the corresponding sampling times before model fitting. Amplitude is constrained to be non-negative, and phase confidence intervals are wrapped to `[0, period)`. When the estimated amplitude is close to zero, phase uncertainty is not well defined and phase standard errors/confidence intervals are returned as `NA`.

## Interpreting results

Most testing functions return a data frame with one row per gene and three columns:

- `pvalue`: raw p-value;
- `BH`: Benjamini-Hochberg adjusted p-value;
- `qvalue`: Storey q-value, or `NA` if the **qvalue** package is unavailable or estimation fails.

A typical workflow uses a false discovery rate threshold such as 0.05:

```r
fdr_level <- 0.05
sig_genes <- rownames(TR_res)[TR_res$BH < fdr_level]
length(sig_genes)
```

For `TDA()` and `TDP()`, tests are applied only to genes that are rhythmic in all conditions according to nominal `TR()` p-values below 0.05. Genes that do not meet this criterion are assigned `NA`.

## Citation

If you use CARhy, please cite:

Huang, W., Menet, J., and Sinha, S. (2026). *CARhy: Comprehensive Analyses of Circadian Rhythms in Transcriptomic Experiments with Multiple Conditions.*

## License

This package is released under the MIT license. See `LICENSE` for details.
