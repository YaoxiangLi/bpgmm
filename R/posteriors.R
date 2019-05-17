updatePostZ = function(m, n, thetaYList){

  tao    = thetaYList@tao
  psy    = thetaYList@psy
  M      = thetaYList@M
  lambda = thetaYList@lambda


  pMat = matrix(NA, m, n)
  ## evaluate density
  dMat = matrix(NA, m, n)

  for(k in 1:m) {
    for(i in 1:n) {
      dMat[k,i] = dmvnorm(X[,i], mean = M[[k]], sigma = psy[[k]] +  lambda[[k]]%*%t(lambda[[k]])
                          ,log = T)
    }
  }

  for(k in 1:m){
    dMat[k,] = dMat[k,] + log(tao[k])
  }

  for(i in 1:n){
    for(k in 1:m){
      pMat[k,i] = calculateRatio(dMat[k,i], dMat[,i])
    }
  }

  ZOneDim = c()
  for(i in 1:n){
    tempProb = as.numeric(pMat[,i])
    ZOneDim[i] = sample(x = 1:m, size = 1, prob = tempProb)
  }
  ZOneDim
}




updatePostThetaY = function(m, n,hparam, thetaYList, ZOneDim, qVec, constraint){

  alpha1 = hparam@alpha1
  alpha2 = hparam@alpha2
  bbeta  = hparam@bbeta
  lambda = thetaYList@lambda
  Y      = thetaYList@Y
  M      = thetaYList@M
  psy    = thetaYList@psy

  ## post for Theta = {tao, M, Lambda, psy}
  CxyList = CalculateCxy(m, n, hparam,thetaYList, ZOneDim, qVec)
  Cxxk = CxyList$Cxxk; Cxyk = CxyList$Cxyk; Cyyk = CxyList$Cyyk;
  Cytytk = CxyList$Cytytk; Cxtytk = CxyList$Cxtytk; CxL1k = CxyList$CxL1k;
  Cxmyk = CxyList$Cxmyk; sumCxmyk = CxyList$sumCxmyk; sumCyyk = CxyList$sumCyyk
  A = CxyList$A; nVec = CxyList$nVec

  Zmat = getZmat(ZOneDim,m,n)


  # post tao
  tao = rdirichlet(1, nVec + ggamma)

  #  post mu
  M = list()
  for(k in 1:m){
    M[[k]] = rmvnorm(1, mean = CxL1k[[k]]*(sum(Zmat[k,]) + alpha1)^(-1) , sigma = (sum(Zmat[k,]) + alpha1)^(-1)* psy[[k]])
  }

  ## lambda; psy
  lambdaPsy = CalculatePostLambdaPsy(alpha1, alpha2, bbeta,CxyList, M, psy, constraint)
  lambda = lambdaPsy$lambda; psy = lambdaPsy$psy

  ## post Y
  D = list()
  for(i in 1:m){
    D[[i]] = t(lambda[[i]]) %*% solve( psy[[i]] +  lambda[[i]]%*%t(lambda[[i]]) )
  }

  Sigma = list()
  for(i in 1:m){
    Sigma[[i]] = diag(qVec[i]) - D[[i]] %*% lambda[[i]]
  }

  Y = list()
  YDvalList = c()
  for(k in 1:m){
    Y[[k]] = matrix(NA, qVec[k], n)
    for(i in 1:n){

      if(Zmat[k,i] == 0){
        Y[[k]][,i] = rmvnorm(1, mean = rep(0,qVec[k]), sigma = diag(qVec[k]))
      }else if(Zmat[k,i] == 1){
        Y[[k]][,i] = rmvnorm(1,mean =  D[[k]]%*%t(X[,i] - M[[k]]), sigma = Sigma[[k]])
      }

    }
  }

  new("ThetaYList", tao=tao, psy=psy, M=M, lambda=lambda, Y=Y)
}


CalculatePostLambdaPsy = function(alpha1, alpha2, bbeta, CxyList, M, psy, constraint){

  Cxxk = CxyList$Cxxk
  Cxyk = CxyList$Cxyk
  Cyyk = CxyList$Cyyk
  Cytytk = CxyList$Cytytk
  Cxtytk = CxyList$Cxtytk
  CxL1k = CxyList$CxL1k
  Cxmyk = CxyList$Cxmyk
  sumCxmyk = CxyList$sumCxmyk
  sumCyyk = CxyList$sumCyyk
  A = CxyList$A
  nVec = CxyList$nVec

  Y      = thetaYList@Y
  M      = thetaYList@M
  psy    = thetaYList@psy
  delta  = hparam@delta
  bbeta  = hparam@bbeta

  ##
  lambda = list()
  if(constraint[1] == T & constraint[2] == T & constraint[3] == T){
    ##model 1
    for(k in 1:m){
      if(k == 1){
        lambda[[k]] = mvtnorm::rmvnorm(1, mean = c(sumCxmyk %*% solve(sumCyyk)),
                              sigma = kronecker(solve(sumCyyk), psy[[k]])
        )
        lambda[[k]] = matrix(lambda[[k]], p, qVec[k])

      }else{
        lambda[[k]] = lambda[[1]]
      }
    }

    # post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k
    tildaLambda = list()
    for(k in 1:m){
      tildaLambda[[k]] = cbind(as.matrix(M[[k]]), lambda[[k]])
    }

    ## post psy
    psy = list()

    shapePara = 0
    ratePara = 0
    for(k in 1:m){
      shapePara = shapePara + p/2 * (nVec[k] + qVec[k] + 2 * delta - 1)
      ratePara = ratePara + 1/2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
                                       + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]])%*%t(tildaLambda[[k]])
                                       + 2 * bbeta * diag(rep(1,p)))

    }
    shapePara = shapePara + 1
    ratePara = sum(ratePara)

    invpsy = rgamma(1, shape = shapePara, rate = ratePara)
    for(k in 1:m){
      psy[[k]] = diag(rep(1/invpsy, p))
    }
    ##model 1 end
  }else if(constraint[1] == T & constraint[2] == T & constraint[3] == F){
    ##model 2
    for(k in 1:m){
      if(k == 1){
        lambda[[k]] = rmvnorm(1, mean = c(sumCxmyk %*% solve(sumCyyk)),
                              sigma = kronecker(solve(sumCyyk), psy[[k]])
        )
        lambda[[k]] = matrix(lambda[[k]], p, qVec[k])

      }else{
        lambda[[k]] = lambda[[1]]
      }
    }

    # post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k
    tildaLambda = list()
    for(k in 1:m){
      tildaLambda[[k]] = cbind(t(M[[k]]), lambda[[k]])
    }

    ## post psy
    psy = list()

    shapePara = 0
    ratePara = 0
    for(k in 1:m){
      shapePara = shapePara + 1/2 * (nVec[k] + qVec[k] + 2 * delta - 1)
      ratePara = ratePara + 1/2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
                                       + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]])%*%t(tildaLambda[[k]] )
                                       + 2 * bbeta * diag(rep(1,p)))

    }
    shapePara = shapePara + 1

    invpsy = c()
    for(j in 1:p){
      invpsy[j] = rgamma(1, shape = shapePara, rate = ratePara[j])
    }
    for(k in 1:m){
      psy[[k]] = diag(1/invpsy)
    }
    ##end model 2
  }else if(constraint[1] == T & constraint[2] == F & constraint[3] == T){
    ##model 3

    sumPhiCxy = 0
    sumPhiCyy = 0

    for(k in 1:m){

      sumPhiCxy = sumPhiCxy + 1/psy[[k]][1,1] * Cxmyk[[k]]
      sumPhiCyy = sumPhiCyy + 1/psy[[k]][1,1] * (Cyyk[[k]] + alpha2 * diag(qVec[k]))
    }

    for(k in 1:m){
      if(k == 1){
        lambda[[k]] =   rmvnorm(1, mean = c(sumPhiCxy %*% solve(sumPhiCyy)),
                                sigma = kronecker(solve(sumPhiCyy), diag(p)))

        lambda[[k]] = matrix(lambda[[k]], p, qVec[k])
      }else{
        lambda[[k]] = lambda[[1]]
      }
    }
    # post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k
    tildaLambda = list()
    for(k in 1:m){
      tildaLambda[[k]] = cbind(t(M[[k]]), lambda[[k]])
    }

    ## post psy
    psy = list()

    shapePara = 0
    ratePara = 0
    for(k in 1:m){
      shapePara = p/2 * (nVec[k] + qVec[k] + 2 * delta - 1) + 1
      ratePara = 1/2 * sum(diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
                                + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]])%*%t(tildaLambda[[k]] )
                                + 2 * bbeta * diag(rep(1,p))))

      invpsy = rgamma(1, shape = shapePara, rate = ratePara)
      psy[[k]] = diag(rep(1/invpsy,p))
    }


    ## end model 3
  }else if(constraint[1] == T & constraint[2] == F & constraint[3] == F){
    ##model 4
    sumVar = 0
    B = 0
    for(k in 1:m){
      sumVar = sumVar + kronecker(Cyyk[[k]] + alpha2 * diag(qVec[k])
                                  , solve(psy[[k]]))
      B = B +  solve(psy[[k]]) %*% Cxmyk[[k]]

    }
    lambdaVar = solve(sumVar)
    lambdaMean = t(c(B)) %*% lambdaVar
    for(k in 1:m){
      if(k == 1){
        lambda[[k]] =   rmvnorm(1, mean = lambdaMean,
                                sigma = lambdaVar)

        lambda[[k]] = matrix(lambda[[k]], p, qVec[k])
      }else{
        lambda[[k]] = lambda[[1]]
      }
    }
    # post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k
    tildaLambda = list()
    for(k in 1:m){
      tildaLambda[[k]] = cbind(t(M[[k]]), lambda[[k]])
    }

    ## post psy
    psy = list()

    shapePara = 0
    ratePara = 0
    for(k in 1:m){
      shapePara = 1/2 * (nVec[k] + qVec[k] + 2 * delta - 1) + 1
      ratePara = 1/2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
                            + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]])%*%t(tildaLambda[[k]] )
                            + 2 * bbeta * diag(rep(1,p)))


      invpsy = c()
      for(j in 1:p){
        invpsy[j] = rgamma(1, shape = shapePara, rate = ratePara[j])
      }
      # invpsy = rgamma(p, shape = shapePara, rate = ratePara)
      psy[[k]] = diag(1/invpsy)
    }

    ##end model 4
  }else if(constraint[1] == F & constraint[2] == T & constraint[3] == T){
    ##model 5
    for(k in 1:m){

      lambda[[k]] = rmvnorm(1, mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]]+  alpha2 * diag(qVec[k]))),
                            sigma = kronecker(solve(sumCyyk), psy[[k]])
      )
      lambda[[k]] = matrix(lambda[[k]], p, qVec[k])


    }

    # post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k
    tildaLambda = list()
    for(k in 1:m){
      tildaLambda[[k]] = cbind(t(M[[k]]), lambda[[k]])
    }

    ## post psy
    psy = list()

    shapePara = 0
    ratePara = 0
    for(k in 1:m){
      shapePara = shapePara + p/2 * (nVec[k] + qVec[k] + 2 * delta - 1)
      ratePara = ratePara + 1/2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
                                       + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]])%*%t(tildaLambda[[k]] )
                                       + 2 * bbeta * diag(rep(1,p)))

    }
    shapePara = shapePara + 1
    ratePara = sum(ratePara)

    invpsy = rgamma(1, shape = shapePara, rate = ratePara)

    for(k in 1:m){
      psy[[k]] = diag(rep(1/invpsy),p)
    }
    ##end model 5
  }else if(constraint[1] == F & constraint[2] == T & constraint[3] == F){
    ##model 6
    for(k in 1:m){

      lambda[[k]] = rmvnorm(1, mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]]+  alpha2 * diag(qVec[k]))),
                            sigma = kronecker(solve(sumCyyk), psy[[k]])
      )
      lambda[[k]] = matrix(lambda[[k]], p, qVec[k])


    }

    # post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k
    tildaLambda = list()
    for(k in 1:m){
      tildaLambda[[k]] = cbind(t(M[[k]]), lambda[[k]])
    }

    ## post psy
    psy = list()

    shapePara = 0
    ratePara = 0
    for(k in 1:m){
      shapePara = shapePara + 1/2 * (nVec[k] + qVec[k] + 2 * delta - 1)
      ratePara = ratePara + 1/2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
                                       + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]])%*%t(tildaLambda[[k]] )
                                       + 2 * bbeta * diag(rep(1,p)))

    }
    shapePara = shapePara + 1

    invpsy = c()
    for(j in 1:p){
      invpsy[j] = rgamma(1, shape = shapePara, rate = ratePara[j])
    }
    for(k in 1:m){
      psy[[k]] = diag(1/invpsy)
    }
    ## end model 6
  }else if(constraint[1] == F & constraint[2] == F & constraint[3] == T){
    ##model 7
    for(k in 1:m){

      lambda[[k]] = rmvnorm(1, mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]]+  alpha2 * diag(qVec[k]))),
                            sigma = kronecker(solve(sumCyyk), psy[[k]])
      )
      lambda[[k]] = matrix(lambda[[k]], p, qVec[k])


    }

    # post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k
    tildaLambda = list()
    for(k in 1:m){
      tildaLambda[[k]] = cbind(t(M[[k]]), lambda[[k]])
    }

    ## post psy
    psy = list()

    shapePara = 0
    ratePara = 0
    for(k in 1:m){
      shapePara = p/2 * (nVec[k] + qVec[k] + 2 * delta - 1) + 1
      ratePara = 1/2 * sum(diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
                                + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]])%*%t(tildaLambda[[k]] )
                                + 2 * bbeta * diag(rep(1,p))))

      invpsy = rgamma(1, shape = shapePara, rate = ratePara)
      psy[[k]] = diag(rep(1/invpsy,p))
    }
    ## end model 7
  }else if(constraint[1] == F & constraint[2] == F & constraint[3] == F){
    ##model 8
    for(k in 1:m){

      lambda[[k]] = rmvnorm(1, mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]] +  alpha2 * diag(qVec[k]))),
                            sigma = kronecker(solve(sumCyyk), psy[[k]])
      )
      lambda[[k]] = matrix(lambda[[k]], p, qVec[k])


    }

    # post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k
    tildaLambda = list()
    for(k in 1:m){
      tildaLambda[[k]] = cbind(t(M[[k]]), lambda[[k]])
    }

    ## post psy
    psy = list()

    shapePara = 0
    ratePara = 0
    for(k in 1:m){
      shapePara = 1/2 * (nVec[k] + qVec[k] + 2 * delta - 1) + 1
      ratePara = 1/2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
                            + tildaLambda[[k]]%*%(Cytytk[[k]] + A[[k]])%*%t(tildaLambda[[k]] )
                            + 2 * bbeta * diag(rep(1,p)))

      invpsy = rgamma(p, shape = shapePara, rate = ratePara)
      psy[[k]] = diag(1/invpsy)
    }

    ##end model 8
  }
  return(list(lambda = lambda,psy = psy) )
}

