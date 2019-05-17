# source("priors.r")
# source("posteriors.r")
# source("pgmmAssist.r")
# source("plotFun.r")
# source("proposalLambda.r")
# source("proposalPsi.r")
# source("jointDensity.r")
# source("HparamClass.R")
# source("paramClass.R")

#' @name bpgmm
#' @aliases bpgmm
#' @title Model-Based Clustering Using Baysian PGMM
#' @description Carries out model-based clustering using parsimonious Gaussian mixture models. MCMC are used for parameter estimation. The RJMCMC is used for model selection.
#' @usage
#'
#'
#'
#'
#'
#'
#'
#' @import stats
#' @param niter the number of iterations.
#' @param burn the number of burn in iterations.
#' @param X the observation matrix with size p * m.
#' @param n the number of observations.
#' @param m the number of clusters.
#' @param delta,ggamma  scaler hyperparameters.
#' @param qVec the vector of the number of factors in each clusters.
#' @param qnew the number of factor for a new cluster.
#' @param constraint the pgmm constraint, a vector of length three with binary entry. For example, c(1,1,1) means the fully constraint model.
#' @param dVec a vector of hyperparameters with length three, shape parameters for alpha1, alpha2 and bbeta respectively.
#' @param sVec a vector of hyperparameters with length three, rate parameters for alpha1, alpha2 and bbeta respectively.

parsimoniousGaussianMixtureModel = function(niter,
                                            burn,
                                            X,
                                            n,
                                            p,
                                            delta,
                                            ggamma,
                                            m,
                                            qVec,
                                            qnew,
                                            constraint,
                                            dVec,
                                            sVec){


  alpha1 = rgamma(1, dVec[1], sVec[1])
  alpha2 = rgamma(1, dVec[2], sVec[2])
  bbeta = rgamma(1, dVec[3], sVec[3])

  hparam = new("Hparam", alpha1=alpha1, alpha2=alpha2, bbeta=bbeta, delta=delta, ggamma=ggamma)

  muBar = apply(X, MARGIN = 1, FUN = mean)

  ## priors
  ZOneDim = kmeans(x = t(X), centers = m)$cluster
  priorList = generatePriorThetaY(m, n, p, muBar, hparam, qVec, ZOneDim, constraint)

  thetaYList = new("ThetaYList",tao = priorList$tao,
                   psy = priorList$psy, M = priorList$M,
                   lambda = priorList$lambda, Y = priorList$Y
                   )

  for(i in 1:burn){
    ZOneDim = updatePostZ(m, n, thetaYList)
    thetaYList = updatePostThetaY(m = m, n = n,hparam, thetaYList, ZOneDim=ZOneDim,  qVec = qVec,constraint = constraint)
    hparam = updateHyperparameter(m, p, qnew, hparam,thetaYList, dVec, sVec)
  }

  ##
  alpha1Vec = c()
  alpha2Vec = c()
  bbetaVec = c()
  taoList = list()
  psyList = list()
  MList = list()
  lambdaList = list()
  YList = list()
  ZmatList = list()
  constraintList = list()
  ##
  diffVec = c()
  lambdaCount = 0

  ## posteriors
  for(h in 1:niter){
    cat("h = ", h, "  =========>\n")
    steps = c("stay", "lambda", "psi1", "psi2")
    currentStep = sample(size = 1, x = steps)
    if(currentStep == "stay"){
      print("stay")
      ZOneDim = updatePostZ(m, n, thetaYList)
      thetaYList = updatePostThetaY(m = m, n = n,hparam, thetaYList, ZOneDim=ZOneDim, qVec = qVec, constraint = constraint)
      hparam = updateHyperparameter(m, p, qnew, hparam, thetaYList, dVec, sVec)

    }else if(currentStep == "lambda"){

      print("lambda")
      proposeConstraint = constraint
      proposeConstraint[1] = (constraint[1] + 1) %% 2

      CxyList = CalculateCxy(m, n, hparam, thetaYList,ZOneDim, qVec)
      newLambda = CalculateProposalLambda(hparam, thetaYList,CxyList, proposeConstraint)

      newthetaYList = thetaYList
      newthetaYList@lambda = newLambda

      newCxyList = CalculateCxy(m, n, hparam,newthetaYList, ZOneDim, qVec)

      oldDensity = likelihood(thetaYList, ZOneDim,qVec,muBar) + evaluatePrior(m, p, muBar,hparam, thetaYList, ZOneDim, qVec,constraint)
      oldLambdaEval = EvaluateProposalLambda(hparam, newthetaYList,newCxyList, constraint, thetaYList@lambda)
      newLambdaEval = EvaluateProposalLambda(hparam, thetaYList,CxyList, proposeConstraint, newthetaYList@lambda)

      ## Gibbs
      newhparam = hparam
      for(i in 1:10){
        newZOneDim = updatePostZ(m, n, newthetaYList)
        newthetaYList = updatePostThetaY(m = m, n = n,newhparam, newthetaYList, ZOneDim=ZOneDim, qVec = qVec, constraint = proposeConstraint)
        newhparam = updateHyperparameter(m, p, qnew, newhparam,newthetaYList, dVec, sVec)
      }

      newDensity = likelihood(newthetaYList, newZOneDim,qVec,muBar) + evaluatePrior(m, p, muBar, newhparam,newthetaYList, newZOneDim, qVec, proposeConstraint)

      numer  = newDensity + oldLambdaEval
      denom = oldDensity + newLambdaEval
      diffVec[h] = numer - denom
      cat("accepting prob = exp ", numer - denom, " \n")

      probAlpha = calculateRatio(numer, denom)
      acceptP = min(1,probAlpha)
      res = rbinom(1, size = 1, prob = acceptP)

      if(res == 1){
        lambdaCount = lambdaCount + 1
        print("lambda success=====>")
        cat("lambda prob = ", probAlpha, "====>\n")

        constraint = proposeConstraint
        thetaYList = newthetaYList
        ZOneDim = newZOneDim
        hparam = newhparam
      }
    }else if(currentStep == "psi1" | currentStep == "psi2"){

      print(currentStep)
      proposeConstraint = constraint
      if(currentStep == "psi1"){
        proposeConstraint[2] = (constraint[2] + 1) %% 2
      }else{
        proposeConstraint[3] = (constraint[3] + 1) %% 2
      }

      #####
      CxyList = CalculateCxy(m, n, hparam, thetaYList , ZOneDim, qVec)
      newPsy = CalculateProposalPsy(hparam, thetaYList, CxyList, proposeConstraint)

      newthetaYList = thetaYList
      newthetaYList@psy = newPsy
      newCxyList = CalculateCxy(m, n, hparam, newthetaYList, ZOneDim, qVec)



      oldDensity = likelihood(thetaYList, ZOneDim,qVec,muBar) + evaluatePrior(m, p, muBar, hparam, thetaYList, ZOneDim, qVec, constraint)
      newDensity = likelihood(newthetaYList, ZOneDim,qVec,muBar) + evaluatePrior(m, p, muBar,hparam, newthetaYList, ZOneDim, qVec, proposeConstraint)

      oldPsyEval = EvaluateProposalPsy(hparam,newthetaYList, newCxyList, constraint, thetaYList@psy)
      newPsyEval = EvaluateProposalPsy(hparam,thetaYList, CxyList, proposeConstraint, newthetaYList@psy)

      numer  = newDensity + oldPsyEval
      denom = oldDensity + newPsyEval
      print(numer - denom)
      diffVec[h] = numer - denom

      probAlpha = calculateRatio(numer, denom)
      acceptP = min(1,probAlpha)
      res = rbinom(1, size = 1, prob = acceptP)

    if(res == 1){
      print("psi success=====>")
      cat("psi prob = ", probAlpha, "====>\n")
      constraint = proposeConstraint
      psy = newPsy
    }
  }

    print(constraint)
    print(ZOneDim)
    ##save
    alpha1Vec[h] = hparam@alpha1
    alpha2Vec[h] = hparam@alpha2
    bbetaVec[h] = hparam@bbeta
    taoList[[h]] = thetaYList@tao
    psyList[[h]] = thetaYList@psy
    MList[[h]] = thetaYList@M
    lambdaList[[h]] = thetaYList@lambda
    YList[[h]] = thetaYList@Y
    ZmatList[[h]] = ZOneDim
    constraintList[[h]] = constraint
  }

  return(list(taoList=taoList, psyList=psyList,MList=MList,lambdaList=lambdaList,
              YList=YList,ZmatList=ZmatList, constraintList = constraintList,lambdaCount = lambdaCount,
              alpha1Vec = alpha1Vec, alpha2Vec = alpha2Vec, bbetaVec = bbetaVec, diffVec = diffVec))
}
