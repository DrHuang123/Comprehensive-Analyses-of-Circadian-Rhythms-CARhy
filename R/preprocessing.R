#' Preprocess RNA-seq count data before analysis
#'
#' Preprocessing of RNA-seq count data is recommended but not required before
#' using this package. Users may apply any suitable normalization workflow.
#' This function provides one optional workflow: low-expression filtering with
#' \code{edgeR::filterByExpr()}, followed by TMM normalization using \pkg{edgeR}.
#' The \pkg{edgeR} package is available from Bioconductor.
#'
#' @param counts A raw count matrix with genes in rows and samples in columns.
#' @param filter Logical; whether to filter lowly expressed genes using  \code{edgeR::filterByExpr()}.
#' @param group Optional group vector passed to  \code{edgeR::filterByExpr()} and  \code{edgeR::DGEList()}.
#' @param design Optional design matrix passed to  \code{edgeR::filterByExpr()}.
#' @param method Normalization method for  \code{edgeR::calcNormFactors()}.
#' @param output One of "logCPM", "CPM", or "DGEList".
#' @param prior.count Prior count used when \code{output = "logCPM"}.
#'
#' @return
#' A normalized expression matrix or an edgeR DGEList object.
#'
#' @references
#' Robinson, M. D., McCarthy, D. J., and Smyth, G. K. (2010).
#' edgeR: a Bioconductor package for differential expression analysis of
#' digital gene expression data. \emph{Bioinformatics}, 26(1), 139-140.
#'
#' @examples
#' # Simulate an RNA-seq count dataset with expression data for 200 genes
#' # under two experimental conditions, A and B, with five replicates per condition.
#'
#' set.seed(123)
#' counts <- matrix(rnbinom(2000, mu = 50, size = 1), nrow = 200, ncol = 10)
#' group <- rep(c("A", "B"), each = 5)
#'
#' expr <- preprocess_counts(
#'  counts,
#'  filter = TRUE,
#'  group = group,
#'  method = "TMM",
#'  output = "logCPM"
#' )
#'
#' head(expr)
#'
#' @export

preprocess_counts <- function(counts,
                              filter = TRUE,
                              group = NULL,
                              design = NULL,
                              method = "TMM",
                              output = c("logCPM", "CPM", "DGEList"),
                              prior.count = 1) {
  output <- match.arg(output)
  counts <- as.matrix(counts)

  if (!is.numeric(counts)) {
    stop("`counts` must be a numeric matrix.")
  }

  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required. Please install it first.")
  }

  y <- edgeR::DGEList(counts = counts, group = group)

  if (filter) {
    if (!is.null(design)) {
      keep <- edgeR::filterByExpr(y, design = design)
    } else {
      keep <- edgeR::filterByExpr(y, group = group)
    }
    y <- y[keep, , keep.lib.sizes = FALSE]
  }

  y <- edgeR::calcNormFactors(y, method = method)

  if (output == "DGEList") {
    return(y)
  }

  if (output == "CPM") {
    return(edgeR::cpm(y, log = FALSE, normalized.lib.sizes = TRUE))
  }

  edgeR::cpm(
    y,
    log = TRUE,
    normalized.lib.sizes = TRUE,
    prior.count = prior.count
  )
}
