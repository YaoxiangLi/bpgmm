#' Bayesian Model-Based Clustering with Parsimonious Gaussian Mixture Models
#'
#' Carries out model-based clustering using parsimonious Gaussian mixture
#' models. MCMC is used for parameter estimation and RJMCMC is used for model
#' selection.
#'
#' The `constraint` argument follows the three-letter PGMM model notation used
#' in Lu, Li, and Love (2021). The first entry indicates whether loading
#' matrices are shared across clusters, the second whether noise covariance
#' matrices are shared across clusters, and the third whether the noise
#' covariance is isotropic within each cluster. Use [model_to_constraint()] to
#' convert model names such as `CCC`, `CCU`, `CUC`, `CUU`, `UCC`, `UCU`, `UUC`,
#' and `UUU` into the numeric vector used internally.
#'
#' @import stats MASS mcmcse pgmm label.switching fabMix
#' @param X the observation matrix with size p * m
#' @param mInit the number of initial clusters
#' @param mVec the range of the number of clusters
#' @param qnew the number of latent factors for a new cluster
#' @param delta scalar hyperparameter for the noise covariance prior
#' @param ggamma scalar hyperparameter used in covariance-structure proposals
#' @param burn the number of burn-in iterations
#' @param niter the number of posterior sampling iterations
#' @param constraint initial PGMM covariance constraint. Use a three-letter
#'   model label such as `"CCC"` or `"UUU"`, or a numeric vector of length
#'   three with binary entries. For example, `c(1, 1, 1)` is `CCC`, the fully
#'   constrained model, and `c(0, 0, 0)` is `UUU`, the fully unconstrained
#'   model.
#' @param dVec a vector of hyperparameters with length three, shape parameters
#'   for alpha1, alpha2 and bbeta respectively
#' @param sVec a vector of hyperparameters with length three, rate parameters
#'   for alpha1, alpha2 and bbeta respectively
#' @param Mstep indicator for RJMCMC model selection on the number of clusters
#' @param Vstep indicator for RJMCMC model selection on covariance structures
#' @param SCind indicator for using split/combine moves in the cluster-number
#'   RJMCMC step
#' @param verbose logical; if `TRUE`, print iteration progress.
#' @name pgmm_rjmcmc
#' @export

pgmm_rjmcmc <- function(X,
                       mInit,
                       mVec,
                       qnew,
                       delta = 2,
                       ggamma = 2,
                       burn = 20,
                       niter = 1000,
                       constraint = c(0, 0, 0),
                       dVec = c(1, 1, 1),
                       sVec = c(1, 1, 1),
                       Mstep = 0,
                       Vstep = 0,
                       SCind = 0,
                       verbose = TRUE) {
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix with variables in rows and observations in columns", call. = FALSE)
  }
  if (any(!is.finite(X))) {
    stop("X must contain only finite numeric values", call. = FALSE)
  }
  if (!isTRUE(length(mVec) == 2L) || any(!is.finite(mVec)) || any(mVec < 1) || mVec[1] > mVec[2]) {
    stop("mVec must be a length-two increasing positive numeric vector", call. = FALSE)
  }
  if (!isTRUE(length(mInit) == 1L) || !is.finite(mInit) || mInit < 1) {
    stop("mInit must be a positive scalar", call. = FALSE)
  }
  if (mInit < mVec[1] || mInit > mVec[2]) {
    stop("mInit must lie within mVec", call. = FALSE)
  }
  if (mInit > ncol(X)) {
    stop("mInit cannot exceed the number of observations in X", call. = FALSE)
  }
  if (!isTRUE(length(qnew) == 1L) || !is.finite(qnew) || qnew < 1) {
    stop("qnew must be a positive scalar", call. = FALSE)
  }
  if (!isTRUE(length(burn) == 1L) || !is.finite(burn) || burn < 0) {
    stop("burn must be a non-negative scalar", call. = FALSE)
  }
  if (!isTRUE(length(niter) == 1L) || !is.finite(niter) || niter < 0) {
    stop("niter must be a non-negative scalar", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE", call. = FALSE)
  }
  if (!isTRUE(length(Mstep) == 1L) || is.na(Mstep) || !(Mstep %in% c(0, 1))) {
    stop("Mstep must be 0 or 1", call. = FALSE)
  }
  if (!isTRUE(length(Vstep) == 1L) || is.na(Vstep) || !(Vstep %in% c(0, 1))) {
    stop("Vstep must be 0 or 1", call. = FALSE)
  }
  if (!isTRUE(length(SCind) == 1L) || is.na(SCind) || !(SCind %in% c(0, 1))) {
    stop("SCind must be 0 or 1", call. = FALSE)
  }

  mInit <- as.integer(mInit)
  mVec <- as.integer(mVec)
  qnew <- as.integer(qnew)
  burn <- as.integer(burn)
  niter <- as.integer(niter)
  Mstep <- as.integer(Mstep)
  Vstep <- as.integer(Vstep)
  SCind <- as.integer(SCind)
  if (is.character(constraint)) {
    constraint <- model_to_constraint(constraint)
  } else {
    constraint <- model_to_constraint(constraint_to_model(constraint))
  }

  n <- ncol(X)
  p <- nrow(X)

  alpha1 <- rgamma(1, dVec[1], sVec[1])
  alpha2 <- rgamma(1, dVec[2], sVec[2])
  bbeta <- rgamma(1, dVec[3], sVec[3])

  hparam <- new("Hparam", alpha1 = alpha1, alpha2 = alpha2, bbeta = bbeta, delta = delta, ggamma = ggamma)

  hparamInit <- hparam

  muBar <- X[, sample.int(n, 1)]

  ## cluster indicator
  clusInd <- rep(0, mVec[2])
  clusInd[1:mInit] <- 1

  ## qinit
  qVec <- rep(0, mVec[2])
  qVec[1:mInit] <- qnew


  ## priors
  ZOneDim <- kmeans(x = t(X), centers = mInit)$cluster
  thetaYList <- generatePriorThetaY(mInit, n, p, muBar, hparam, qVec, ZOneDim, constraint)

  ## burn in
  for (i in seq_len(burn)) {
    MCMCobj <- stayMCMCupdate(X, thetaYList, ZOneDim, hparam, qVec, qnew, dVec, sVec, constraint, clusInd)
    ZOneDim <- MCMCobj$ZOneDim
    thetaYList <- MCMCobj$thetaYList
    hparam <- MCMCobj$hparam
    hparam@alpha2 <- max(0.01, hparam@alpha2)
  }

  thetaYList <- clearCurrentThetaYlist(thetaYList, clusInd, mVec[2])
  ##
  alpha1Vec <- c()
  alpha2Vec <- c()
  bbetaVec <- c()
  taoList <- list()
  psyList <- list()
  MList <- list()
  lambdaList <- list()
  YList <- list()
  ZmatList <- list()
  constraintList <- list()
  clusIndList <- list()
  ##

  for (h in seq_len(niter)) {
    if (verbose) {
      cat("iter = ", h, "======>\n")
    }

    ## choose m or choose v
    if (Mstep == 1) {
      MCMCobj <- MstepRJMCMCupdate(X, muBar, p, thetaYList, ZOneDim, hparam, hparamInit, qVec, qnew, dVec, sVec, constraint, clusInd, mVec, "BD")
      ZOneDim <- MCMCobj$ZOneDim
      thetaYList <- MCMCobj$thetaYList
      hparam <- MCMCobj$hparam
      qVec <- MCMCobj$qVec
      clusInd <- MCMCobj$clusInd
      ##
      if (SCind == 1) {
        MCMCobj <- MstepRJMCMCupdate(X, muBar, p, thetaYList, ZOneDim, hparam, hparamInit, qVec, qnew, dVec, sVec, constraint, clusInd, mVec, "SC")
        ZOneDim <- MCMCobj$ZOneDim
        thetaYList <- MCMCobj$thetaYList
        hparam <- MCMCobj$hparam
        qVec <- MCMCobj$qVec
        clusInd <- MCMCobj$clusInd
      }
    }

    if (Vstep == 1) {
      MCMCobj <- VstepRJMCMCupdate(X, muBar, p, thetaYList, ZOneDim, hparam, hparamInit, qVec, qnew, ggamma, dVec, sVec, constraint, clusInd)
      ZOneDim <- MCMCobj$ZOneDim
      thetaYList <- MCMCobj$thetaYList
      hparam <- MCMCobj$hparam
      constraint <- MCMCobj$constraint
    }

    # stay step
    MCMCobj <- stayMCMCupdate(X, thetaYList, ZOneDim, hparam, qVec, qnew, dVec, sVec, constraint, clusInd)
    ZOneDim <- MCMCobj$ZOneDim
    thetaYList <- MCMCobj$thetaYList
    hparam <- MCMCobj$hparam
    hparam@alpha2 <- max(0.01, hparam@alpha2)

    ## save
    clusIndList[[h]] <- clusInd
    alpha1Vec[h] <- hparam@alpha1
    alpha2Vec[h] <- hparam@alpha2
    bbetaVec[h] <- hparam@bbeta
    taoList[[h]] <- thetaYList@tao
    psyList[[h]] <- thetaYList@psy
    MList[[h]] <- thetaYList@M
    lambdaList[[h]] <- thetaYList@lambda
    YList[[h]] <- thetaYList@Y
    ZmatList[[h]] <- ZOneDim
    constraintList[[h]] <- constraint
  }

  list(
    taoList = taoList, psyList = psyList, MList = MList, lambdaList = lambdaList,
    YList = YList, ZmatList = ZmatList, constraintList = constraintList,
    alpha1Vec = alpha1Vec, alpha2Vec = alpha2Vec, bbetaVec = bbetaVec,
    clusIndList = clusIndList
  )
}

#' @rdname pgmm_rjmcmc
#' @export
pgmmRJMCMC <- function(X,
                       mInit,
                       mVec,
                       qnew,
                       delta = 2,
                       ggamma = 2,
                       burn = 20,
                       niter = 1000,
                       constraint = c(0, 0, 0),
                       dVec = c(1, 1, 1),
                       sVec = c(1, 1, 1),
                       Mstep = 0,
                       Vstep = 0,
                       SCind = 0,
                       verbose = TRUE) {
  .Deprecated("pgmm_rjmcmc")
  pgmm_rjmcmc(
    X = X,
    mInit = mInit,
    mVec = mVec,
    qnew = qnew,
    delta = delta,
    ggamma = ggamma,
    burn = burn,
    niter = niter,
    constraint = constraint,
    dVec = dVec,
    sVec = sVec,
    Mstep = Mstep,
    Vstep = Vstep,
    SCind = SCind,
    verbose = verbose
  )
}
