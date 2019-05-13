# updatePostZ = function(m, n, thetaYList){
#
#   tao = thetaYList@tao
#   psy = thetaYList@psy
#     M = thetaYList@M
#   lambda = thetaYList@lambda
#
#
#   pMat = matrix(NA, m, n)
#   ## evaluate density
#   dMat = matrix(NA, m, n)
#
#   for(k in 1:m){
#     for(i in 1:n){
#       dMat[k,i] = dmvnorm(X[,i], mean = M[[k]], sigma = psy[[k]] +  lambda[[k]]%*%t(lambda[[k]])
#                           ,log = T)
#     }
#   }
#
#   for(k in 1:m){
#     dMat[k,] = dMat[k,] + log(tao[k])
#   }
#
#   for(i in 1:n){
#     for(k in 1:m){
#       pMat[k,i] = calculateRatio(dMat[k,i], dMat[,i])
#     }
#   }
#
#   ZOneDim = c()
#   for(i in 1:n){
#     tempProb = as.numeric(pMat[,i])
#     ZOneDim[i] = sample(x = 1:m, size = 1, prob = tempProb)
#   }
#   ZOneDim
# }
#
#
#' @export
r_CalculateCxy = function(m, n, hparam, thetaYList, ZOneDim, qVec, X){

  alpha1 = hparam@alpha1
  alpha2 = hparam@alpha2
  bbeta = hparam@bbeta

  Y = thetaYList@Y
  lambda = thetaYList@lambda
  M = thetaYList@M
  psy = thetaYList@psy

  Zmat = getZmat(ZOneDim,m,n)

  A = list()
  nVec = c()
  for(k in 1:m){
    nVec[k] = sum(Zmat[k,])
    A[[k]] = diag(c(alpha1, alpha2 * rep(1, qVec[k])))
  }

  Cxxk = list()
  Cxyk = list()
  Cyyk = list()
  Cytytk = list()
  Cxtytk = list()
  CxL1k = list()
  Cxmyk = list()

  for(k in 1:m){
    Cxxk[[k]] = 0
    Cxyk[[k]] = 0
    Cyyk[[k]] = 0
    Cytytk[[k]] = 0
    CxL1k[[k]] = 0
    Cxmyk[[k]] = 0
    Cxtytk[[k]] = 0
    for(i in 1:n){
      Cxxk[[k]] = Cxxk[[k]] + Zmat[k,i] * X[,i] %*% t(X[,i])
      Cxyk[[k]] = Cxyk[[k]] + Zmat[k,i] * X[,i] %*% t(c(Y[[k]][,i]))
      Cyyk[[k]] = Cyyk[[k]] +  Zmat[k,i] * c(Y[[k]][,i]) %*% t(c(Y[[k]][,i]))
      Cxtytk[[k]] = Cxtytk[[k]] + Zmat[k,i] * X[,i] %*% t(c(1,Y[[k]][,i]))
      Cytytk[[k]] = Cytytk[[k]] + Zmat[k,i] * c(1,Y[[k]][,i]) %*% t(c(1,Y[[k]][,i]))

      Cxmyk[[k]] = Cxmyk[[k]] + Zmat[k,i] * c(X[,i] - M[[k]]) %*% t(c(Y[[k]][,i]))
      CxL1k[[k]] = CxL1k[[k]] + Zmat[k,i] * (X[,i] - lambda[[k]]%*%c(Y[[k]][,i]))
    }
  }

  sumCxmyk = 0
  sumCyyk = 0
  for(k in 1:m){
    sumCxmyk = sumCxmyk + Cxmyk[[k]]
    sumCyyk = sumCyyk + Cyyk[[k]] + alpha2*diag(qVec[k])
  }
  return(list(A = A,
              nVec = t(t(nVec)),
              Cxxk = Cxxk,
              Cxyk = Cxyk,
              Cyyk = Cyyk,
              Cytytk = Cytytk,
              Cxtytk = Cxtytk,
              CxL1k = CxL1k,
              Cxmyk = Cxmyk,
              sumCxmyk = sumCxmyk,
              sumCyyk = sumCyyk))
}
#
#
# updatePostThetaY = function(m, n,hparam, thetaYList, ZOneDim, qVec, constraint){
#
#   alpha1 = hparam@alpha1
#   alpha2 = hparam@alpha2
#   bbeta = hparam@bbeta
#   lambda = thetaYList@lambda
#   Y = thetaYList@Y
#   M = thetaYList@M
#   psy = thetaYList@psy
#
#   ## post for Theta = {tao, M, Lambda, psy}
#   CxyList = CalculateCxy(m, n, hparam,thetaYList, ZOneDim, qVec)
#   Cxxk = CxyList$Cxxk; Cxyk = CxyList$Cxyk; Cyyk = CxyList$Cyyk;
#   Cytytk = CxyList$Cytytk; Cxtytk = CxyList$Cxtytk; CxL1k = CxyList$CxL1k;
#   Cxmyk = CxyList$Cxmyk; sumCxmyk = CxyList$sumCxmyk; sumCyyk = CxyList$sumCyyk
#   A = CxyList$A; nVec = CxyList$nVec
#
#   Zmat = getZmat(ZOneDim,m,n)
#
#
#   # post tao
#   tao = rdirichlet(1, nVec + ggamma)
#
#   #  post mu
#   M = list()
#   for(k in 1:m){
#     M[[k]] = rmvnorm(1, mean = CxL1k[[k]]*(sum(Zmat[k,]) + alpha1)^(-1) , sigma = (sum(Zmat[k,]) + alpha1)^(-1)* psy[[k]])
#   }
#
#   ## lambda; psy
#   lambdaPsy = CalculatePostLambdaPsy(alpha1, alpha2, bbeta,CxyList, M, psy, constraint)
#   lambda = lambdaPsy$lambda; psy = lambdaPsy$psy
#
#   ## post Y
#   D = list()
#   for(i in 1:m){
#     D[[i]] = t(lambda[[i]]) %*% solve( psy[[i]] +  lambda[[i]]%*%t(lambda[[i]]) )
#   }
#
#   Sigma = list()
#   for(i in 1:m){
#     Sigma[[i]] = diag(qVec[i]) - D[[i]] %*% lambda[[i]]
#   }
#
#   Y = list()
#   YDvalList = c()
#   for(k in 1:m){
#     Y[[k]] = matrix(NA, qVec[k], n)
#     for(i in 1:n){
#
#       if(Zmat[k,i] == 0){
#         Y[[k]][,i] = rmvnorm(1, mean = rep(0,qVec[k]), sigma = diag(qVec[k]))
#       }else if(Zmat[k,i] == 1){
#         Y[[k]][,i] = rmvnorm(1,mean =  D[[k]]%*%t(X[,i] - M[[k]]), sigma = Sigma[[k]])
#       }
#
#     }
#   }
#
#   new("ThetaYList", tao=tao, psy=psy, M=M, lambda=lambda, Y=Y)
#   # return(list(tao=tao, psy=psy, M=M, lambda=lambda, Y=Y))
# }
#
