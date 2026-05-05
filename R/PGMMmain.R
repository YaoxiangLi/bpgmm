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
#' @param constraint initial PGMM covariance constraint, a numeric vector of
#'   length three with binary entries. For example, `c(1, 1, 1)` is `CCC`, the
#'   fully constrained model, and `c(0, 0, 0)` is `UUU`, the fully unconstrained
#'   model.
#' @param dVec a vector of hyperparameters with length three, shape parameters
#'   for alpha1, alpha2 and bbeta respectively
#' @param sVec a vector of hyperparameters with length three, rate parameters
#'   for alpha1, alpha2 and bbeta respectively
#' @param Mstep indicator for RJMCMC model selection on the number of clusters
#' @param Vstep indicator for RJMCMC model selection on covariance structures
#' @param SCind indicator for using split/combine moves in the cluster-number
#'   RJMCMC step
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
                       SCind = 0) {
  n <- ncol(X)
  p <- nrow(X)

  alpha1 <- rgamma(1, dVec[1], sVec[1])
  alpha2 <- rgamma(1, dVec[2], sVec[2])
  bbeta <- rgamma(1, dVec[3], sVec[3])

  hparam <- new("Hparam", alpha1 = alpha1, alpha2 = alpha2, bbeta = bbeta, delta = delta, ggamma = ggamma)

  hparamInit <- hparam

  muBar <- X[, sample(1:n, 1)]

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
  for (i in 1:burn) {
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

  for (h in 1:niter) {
    cat("iter = ", h, "======>\n")

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
                       SCind = 0) {
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
    SCind = SCind
  )
}
