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
- Reports raw p-values, Benjamini-Hochberg adjusted p-values, and Storey q-values.
- Optional RNA-seq count preprocessing using low-expression filtering and TMM normalization via **edgeR**.

## Installation

```r
install.packages("remotes")
remotes::install_github("DrHuang123/Comprehensive-Analyses-of-Circadian-Rhythms-CARhy")
```

## Optional packages

If you would like to use Storey q-value estimation, install the **qvalue** package:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("qvalue")
```

If you would like to use the optional RNA-seq count preprocessing workflow in `preprocess_counts()`, install **edgeR**:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("edgeR")
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

## Documentation

For more details about each function, see the corresponding help page in R:

```r
?preprocess_counts
?TR
?TDR
?TDA
?TDP
?TDM
?params_output
```

## Citation

If you use CARhy, please cite:

Huang, W., Menet, J., and Sinha, S. (2026). *CARhy: Comprehensive Analyses of Circadian Rhythms in Transcriptomic Experiments with Multiple Conditions.*

## License

This package is released under the MIT license. See `LICENSE` for details.
