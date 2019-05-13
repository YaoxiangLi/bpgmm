#' @title PriorThetaY list
#' @importFrom gtools rdirichlet
#' @importFrom mvtnorm rmvnorm
#' @import stats
#'
#' @export
generatePriorThetaY = function(m,
                               n,
                               p,
                               muBar,
                               hparam,
                               qVec,
                               ZOneDim,
                               constraint) {
  ggamma <- hparam@ggamma
  delta  <- hparam@delta
  bbeta  <- hparam@bbeta
  alpha1 <- hparam@alpha1
  alpha2 <- hparam@alpha2

  # prior tao
  if (m == 1) {
    tao <- rbeta(n = 1, shape1 = 1, shape2 = m)
  } else{
    tao <- gtools::rdirichlet(n = 1, alpha = rep(ggamma, m))
  }

  # prior psy
  psy = generatePriorPsi(p, m, delta, bbeta, constraint)

  # prior M
  M = list()
  for (i in 1:m) {
    M[[i]] = c(mvtnorm::rmvnorm(n =  1, mean = muBar, sigma = 1 / alpha1 * psy[[i]]))
  }

  # prior lambda
  lambda = generatePriorLambda(p, m, alpha2, qVec, psy, constraint)


  Zmat = getZmat(ZOneDim, m, n)

  # post Y
  Y = list()
  for (k in 1:m) {
    Y[[k]] = matrix(NA, qVec[k], n)
    for (i in 1:n) {
      Y[[k]][, i] = mvtnorm::rmvnorm(1, mean = rep(0, qVec[k]), sigma = diag(qVec[k]))
    }
  }

  new("ThetaYList",
      tao = tao,
      psy = psy,
      M = M,
      lambda = lambda,
      Y = Y)

  # return(list(
  #   tao    = tao,
  #   psy    = psy,
  #   M      = M,
  #   lambda = lambda,
  #   Y      = Y
  # ))
}

#' evaluate Prior
#'
#'
#' @export

evaluatePrior = function(m,p, muBar,hparam, thetaYList, ZOneDim,qVec, constraint){
  ## product of  2.2, 2.3, 2.6, 2.7, 2.8, 2.9, all in log scale

  ggamma = hparam@ggamma
  alpha1 = hparam@alpha1
  alpha2 = hparam@alpha2
  bbeta  = hparam@bbeta
  delta  = hparam@delta
  ## 2.2: Y
  Yval = 0
  for(k in 1:m){
    Yval = Yval + sum(dnorm(thetaYList@Y[[k]], log = T))
  }
  ## 2.3: Z
  adjustTao = thetaYList@tao/sum(thetaYList@tao)
  Zval = log(sum(adjustTao[ZOneDim]))
  ## 2.6: tao
  taoVal = log(gtools::ddirichlet(x = adjustTao, alpha = rep(ggamma, m)))

  ## 2.7: M
  Mval = 0
  for(k in 1:m){
    Mval = Mval + mvtnorm::dmvnorm(c(thetaYList@M[[k]]), mean = c(muBar), sigma = 1/alpha1 * thetaYList@psy[[k]], log = T)
  }

  ## 2.8: lambda
  lambdaVal = evaluatePriorLambda(p, m, alpha2, qVec,thetaYList@psy, thetaYList@lambda, constraint)

  ## 2.9: psy
  psyVal = evaluatePriorPsi(thetaYList@psy, p, m, delta, bbeta, constraint)

  totalVal = sum(Yval + Zval + taoVal + Mval + lambdaVal + psyVal)
  return(totalVal)

}

#' @export
#'
#'
#'
generatePriorPsi = function(p, m, delta, bbeta, constraint){
  psy = list()


  if(constraint[2] == T & constraint[3] == T){

    for(i in 1:m){
      if(i == 1){
        psyValue = 1/rgamma(1, shape = delta, rate = bbeta)
        psy[[i]] = diag(rep(psyValue, p))
      }else{
        psy[[i]] = psy[[1]]
      }
    }

  }else if(constraint[2] == T & constraint[3] == F){
    for(i in 1:m){
      if(i == 1){
        psy[[i]] = diag(1/rgamma(p, shape = delta, rate = bbeta))
      }else{
        psy[[i]] =  psy[[1]]
      }
    }
  }else if(constraint[2] == F & constraint[3] == T){

    for(i in 1:m){
      psyValue = 1/rgamma(1, shape = delta, rate = bbeta)
      psy[[i]] = diag(rep(psyValue, p))
    }

  }else{
    for(i in 1:m){
      psy[[i]] = diag(1/rgamma(p, shape = delta, rate = bbeta))
    }
  }
  return(psy)
}

#'
#'
#'
#' @export
evaluatePriorPsi = function(psy, p, m, delta, bbeta, constraint){
  psyeval = 0


  if(constraint[2] == T & constraint[3] == T){

    for(i in 1:m){
      if(i == 1){
        psyValue = 1/psy[[i]][1,1]
        psyeval = psyeval + dgamma(psyValue, shape = delta, rate = bbeta, log = T)
      }
    }

  }else if(constraint[2] == T & constraint[3] == F){
    for(i in 1:m){
      if(i == 1){
        psyValue = 1/diag(psy[[i]])
        psyeval = psyeval + sum(dgamma(psyValue, shape = delta, rate = bbeta, log = T))
      }
    }
  }else if(constraint[2] == F & constraint[3] == T){

    for(i in 1:m){
      psyValue = 1/psy[[i]][1,1]
      psyeval = psyeval + dgamma(psyValue, shape = delta, rate = bbeta, log = T)
    }

  }else{
    for(i in 1:m){
      psyValue = 1/diag(psy[[i]])
      psyeval = psyeval + sum(dgamma(psyValue, shape = delta, rate = bbeta, log = T))
    }
  }
  return(psyeval)
}

#' @export
generatePriorLambda = function(p, m, alpha2, qVec,psy, constraint){

  lambda = list()
  if(constraint[1] == T){
    for(k in 1:m){
      if(k == 1){
        qk = qVec[k]
        lambda[[k]] = matrix(0, p, qk)
        for(j in 1:qk){
          lambda[[k]][,j] = mvtnorm::rmvnorm(1, rep(0,p), 1/alpha2 * psy[[k]])
        }
      }else{
        lambda[[k]] = lambda[[1]]
      }
    }
  }else{
    for(k in 1:m){
      qk = qVec[k]
      lambda[[k]] = matrix(0, p, qk)
      for(j in 1:qk){
        lambda[[k]][,j] = mvtnorm::rmvnorm(1, rep(0,p), 1/alpha2 * psy[[k]])
      }
    }
  }
  return(lambda)
}


#' @export
evaluatePriorLambda = function(p, m, alpha2, qVec,psy, lambda, constraint){

  evallambda = 0
  if(constraint[1] == T){
    for(k in 1:m){
      if(k == 1){
        # qk = qVec[k]
        for(j in 1:qk){
          evallambda = evallambda + mvtnorm::dmvnorm(lambda[[k]][,j], rep(0,p), 1/alpha2 * psy[[k]], log = T)
        }
      }
    }
  }else{
    for(k in 1:m){
      qk = qVec[k]
      for(j in 1:qk){
        evallambda = evallambda + mvtnorm::dmvnorm(lambda[[k]][,j], rep(0,p), 1/alpha2 * psy[[k]], log = T)
      }
    }
  }
  return(evallambda)
}
