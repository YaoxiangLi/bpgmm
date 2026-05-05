#' Summarize RJMCMC Samples from a Bayesian PGMM Fit
#'
#' Summarizes posterior samples from [pgmmRJMCMC()] into the modal allocation,
#' posterior counts for the number of clusters, posterior counts for the eight
#' PGMM covariance-constraint models, and optionally the adjusted Rand index
#' against a known reference partition.
#'
#' @param pgmmResList Result list from [pgmmRJMCMC()].
#' @param trueCluster Optional true or reference cluster allocation.
#'
#' @return A list with `Zalloc`, `nCluster`, `nConstraint`, and optionally
#'   `ari`.
#' @importFrom mclust adjustedRandIndex
#' @export
summarizePgmmRJMCMC <- function(pgmmResList, trueCluster = NULL) {
  Zalloc <- sumerizeZ(pgmmResList$ZmatList)

  nCluster <- table(sapply(pgmmResList$ZmatList, function(x) {
    length(unique(x))
  }))

  nConstraint <- pgmmResList$constraintList
  nConstraint <- listToStrVec(nConstraint)
  nConstraint <- table(nConstraint, dnn = "")



  sumRes <- list(Zalloc = Zalloc, nCluster = nCluster, nConstraint = nConstraint)

  if (!is.null(trueCluster)) {
    ari <- adjustedRandIndex(trueCluster, Zalloc)
    sumRes$ari <- ari
  }

  sumRes
}

#' @rdname summarizePgmmRJMCMC
#' @export
summerizePgmmRJMCMC <- function(pgmmResList, trueCluster = NULL) {
  summarizePgmmRJMCMC(pgmmResList, trueCluster)
}
