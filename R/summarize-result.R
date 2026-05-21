#' Summarize RJMCMC Samples from a Bayesian PGMM Fit
#'
#' Summarizes posterior samples from [pgmm_rjmcmc()] into the modal allocation,
#' posterior counts for the number of clusters, posterior counts for the eight
#' PGMM covariance-constraint models, and optionally the adjusted Rand index
#' against a known reference partition.
#'
#' @param pgmmResList Result list from [pgmm_rjmcmc()].
#' @param trueCluster Optional true or reference cluster allocation.
#'
#' @return A list with `Zalloc`, `nCluster`, `nConstraint`, and optionally
#'   `ari`.
#' @importFrom mclust adjustedRandIndex
#' @name summarize_pgmm_rjmcmc
#' @export
summarize_pgmm_rjmcmc <- function(pgmmResList, trueCluster = NULL) {
  if (!is.list(pgmmResList)) {
    stop("pgmmResList must be a result list from pgmm_rjmcmc()", call. = FALSE)
  }
  if (!is.list(pgmmResList$ZmatList) || length(pgmmResList$ZmatList) == 0L) {
    stop("pgmmResList$ZmatList must contain at least one allocation sample", call. = FALSE)
  }
  if (!is.list(pgmmResList$constraintList) || length(pgmmResList$constraintList) != length(pgmmResList$ZmatList)) {
    stop("pgmmResList$constraintList must match pgmmResList$ZmatList", call. = FALSE)
  }

  Zalloc <- summarize_allocations(pgmmResList$ZmatList)

  nCluster <- table(sapply(pgmmResList$ZmatList, function(x) {
    length(unique(x))
  }))

  nConstraint <- pgmmResList$constraintList
  nConstraint <- constraint_list_to_models(nConstraint)
  nConstraint <- table(nConstraint, dnn = "")



  sumRes <- list(Zalloc = Zalloc, nCluster = nCluster, nConstraint = nConstraint)

  if (!is.null(trueCluster)) {
    ari <- adjustedRandIndex(trueCluster, Zalloc)
    sumRes$ari <- ari
  }

  sumRes
}

#' @rdname summarize_pgmm_rjmcmc
#' @export
summarizePgmmRJMCMC <- function(pgmmResList, trueCluster = NULL) {
  .Deprecated("summarize_pgmm_rjmcmc")
  summarize_pgmm_rjmcmc(pgmmResList, trueCluster)
}

#' @rdname summarize_pgmm_rjmcmc
#' @export
summerizePgmmRJMCMC <- function(pgmmResList, trueCluster = NULL) {
  .Deprecated("summarize_pgmm_rjmcmc")
  summarize_pgmm_rjmcmc(pgmmResList, trueCluster)
}
