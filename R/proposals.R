#
#============== Set of functions for full conditionals and proposal distributions
#


sigmaFc <- function(par, mu, nu, tau)
  # Full conditional for sigmaAlpha and sigmaBeta
{
  K <- length(par)
  sigmaNew <- 1/rgamma( 1, shape = (nu + K + 1)/2,
                        scale = ( tau * ( 1 + sum((par - mu)^2) ) + mu*mu ) / (2*tau) )
  return(sigmaNew)
}


muFc <- function(par, sigma, tau, mu, lBound, uBound)
  # Full conditional for muAlpha and muBeta
{
  K <- length(par)
  den <- 1 + tau*K
  muNew <- RcppTN::rtn(.mean = (tau * sum(par) + mu) / den,
                       .sd = sqrt( (tau*sigma)/den ),
                       .low = lBound, .high = uBound)
  return(muNew)
}


alphaBetaProp_k <- function(k, Y, n, arcSumY, alpha, beta, gammaTheta, distZ,
                            muAlpha, sigmaAlpha, muBeta, sigmaBeta,
                            boundA, lambdaCovSum, ind, noSendRec)
  # Proposal for intercept and slope parameter for each view
{
  # alpha -----------------------
  arg <- muAlpha*gammaTheta - beta*distZ - lambdaCovSum   # if no covariates, lambdaCovSum = 0
  Ek <- if ( noSendRec ) arcSumY else sum( (Y*gammaTheta)[ind] )
  expArg <- exp(arg)
  t1 <- expArg / ( (1 + expArg)^2 ) * ( gammaTheta^2 )
  t2 <- expArg / (1 + expArg) * gammaTheta
  diag(t1) <- diag(t2) <- 0
  varAlphaNew_k <- 1/( sum(t1[ind]) + 1/sigmaAlpha )
  meanAlphaNew_k <- varAlphaNew_k * ( Ek - sum(t2[ind]) ) + muAlpha
  varAlphaNew_k <- varAlphaNew_k + n/200
  alphaNew_k <-
    RcppTN::rtn(.mean = meanAlphaNew_k, .sd = sqrt(varAlphaNew_k),
                .low = boundA, .checks = FALSE)

  # beta ------------------------
  arg <- alpha*gammaTheta - muBeta*distZ - lambdaCovSum
  expArg <- exp(arg)
  tmp <- distZ^2 * (expArg/(1+expArg)^2)
  diag(tmp) <- 0
  varBetaNew_k <- 1 / ( sum(tmp[ind]) + (1/sigmaBeta) )
  tmp <- distZ*( ( expArg/(1+expArg) ) - Y )
  diag(tmp) <- 0
  meanBetaNew_k <- varBetaNew_k*sum(tmp[ind]) + muBeta
  varBetaNew_k <- varBetaNew_k + n/200
  betaNew_k <-
    RcppTN::rtn(.mean = meanBetaNew_k, .sd = sqrt(varBetaNew_k),
                .low = 0, .checks = FALSE)
  return( matrix( c(alphaNew_k, meanAlphaNew_k, varAlphaNew_k,
                    betaNew_k, meanBetaNew_k, varBetaNew_k), 1, 6 ) )
  # returns a matrix whose columns are indexed by k
}


zProp_i <- function(i, Y_i, n, K, D, z, distZ_i, alpha, beta, ind,
                    gammaTheta_i, lambdaCovSum_i)
  # Proposal for latent positions for each node
{
  const1 <- matrix(alpha/beta, n, K, byrow = TRUE)*gammaTheta_i
  const2 <- if ( length(lambdaCovSum_i) > 1 ) {
    matrix(lambdaCovSum_i, n, K) / matrix(beta, n, K, byrow = TRUE)
  } else 0
  W <- ( (const1 - const2) > matrix(distZ_i, n, K) )*1    # ausiliary variable
  dYW <- (Y_i - W) #dim:n K

  varZ_i <- 1/(2*sum( sapply(1:K, function(x) sum(abs(dYW[,x])[ind[,x]] ) )* beta ) + 1)    # D dimensional matrix

  tmp <-vapply(1:K, function(k) dYW[,k]*z*beta[k] , numeric(D*n)) ##dim (D*n) K
  tmp <- vapply( 1:K, function(k)  sapply( seq(1,D*n, by =n),
                                           function(x)  sum( tmp[x:(x+n-1),k][ind[,k]])), numeric(D)) ##dim D *K

  meanZ_i <- if( D > 1) rowSums( tmp*2*varZ_i) else sum( tmp*2*varZ_i)
  varZ_i <- varZ_i + n/100   # inflate variance
  zNew_i <- MASS::mvrnorm(1, meanZ_i, diag(varZ_i, D))

  return( list(zNew_i = zNew_i, meanZ_i = meanZ_i, varZ_i = varZ_i) )
}


lambdaProp_l <- function(l, Y, n, K, nC, lambdaCov, lambdaCovSum, muLambda_l, sigmaLambda_l,
                         alpha, beta, noSendRec, gammaTheta, distZ, ind)
  # Proposal for (lth) coefficient parameter for covariates
{
  out <- vapply(1:K, function(k) {
    gammaTheta_k <- if ( noSendRec ) gammaTheta else gammaTheta[,,k]
    arg <- exp( alpha[k]*gammaTheta_k - beta[k]*distZ - (lambdaCovSum - lambdaCov[,,l]) )
    tmp <- (lambdaCov[,,l])^2 * arg/(1 + arg)^2
    diag(tmp) <- 0
    somma1 <- sum(tmp[ind[,,k]])
    tmp <- (arg/(1 + arg) - Y[,,k]) * lambdaCov[,,l]
    diag(tmp) <- 0
    somma2 <- sum(tmp[ind[,,k]])
    return( cbind(somma1, somma2) )
  }, numeric(2) )

  varLambdaNew_l <- 1/( sum(out[1,]) + 1/sigmaLambda_l[l] )
  meanLambdaNew_l <- sum(out[2,])*varLambdaNew_l + muLambda_l[l]

  lambdaNew_l <- RcppTN::rtn(.mean = meanLambdaNew_l, .sd = sqrt(varLambdaNew_l),
                             .low = 0, .checks = FALSE)
  return( matrix(c(lambdaNew_l, meanLambdaNew_l, varLambdaNew_l), 1,3) )
  # returns a matrix whose columns are indexed by l (prop value, mean and variance for each l)
}


gammaThetaProp_i <- function(par_i_old, Y_i, n, K, distZ_i,i,
                             alpha, beta, gt, gammaTheta_i, otherPar,
                             lambdaCovSum_i, effect, effect1,         # effect1 is for the refernce nodes
                             meanParOld_i = 0, varParOld_i = 0, ind)  # effect is either N or C or V
{
  if( is.null(effect)){ parProp = par_i_old
  accRatio = 0
  muProp = varProp = NULL}else{
    if( effect == "const"){
      if(effect1[1] == i){    # fixed effects
        parProp = par_i_old
        accRatio = 0
        muProp =  meanParOld_i
        varProp = varParOld_i
      }else{
        arg <- exp(sapply(1:K, function(x)
          alpha[x] *gammaTheta_i[,x] - beta[x]*distZ_i -lambdaCovSum_i)) # dimension is n x K

        tmp <- arg/((1 +arg)^2)
        tmp <- sapply(1:K, function(x) sum(tmp[,x][ind[,x]] ))
        varProp <- rep(( (1/(gt^2)) *sum( (alpha^2) *tmp))^(-1), K) # 1*K

        tmp <- Y_i -arg/(1 +arg)
        tmp <- sapply(1:K, function(x) sum(tmp[,x][ind[,x]] ))
        muProp <- rep(sum(alpha/gt * tmp)*varProp[1] +(par_i_old[1]), K) # 1*K

        parProp <- rep( RcppTN::rtn(.mean = muProp[1], .sd = sqrt(varProp[1]), # 1
                                    .low = -1, .high = 1, .checks = FALSE), K)
        accRatio = log(RcppTN::dtn(par_i_old[1], .mean = muProp[1], .sd = sqrt(varProp[1]),
                                   .low = -1, .high = 1,.checks = FALSE)) -
          log(RcppTN::dtn(parProp[1], .mean = meanParOld_i[1], .sd = sqrt(varParOld_i[1]),
                          .low = -1, .high = 1, .checks = FALSE))
      }
    } # close const
    if( effect == "var"){
      nullEffect <- which(is.na(par_i_old) )
      fixedEff <- which(effect1 == i)
      fix_nullEff <- sort(c(nullEffect, fixedEff ))
      freeEff <- setdiff(1:K, fix_nullEff)
      lfE <- length(freeEff)
      arg <- exp(sapply(1:K, function(x)
        alpha[x] *gammaTheta_i[,x] -beta[x]*distZ_i -lambdaCovSum_i)) # dimension is n x K
      if( lfE == 0){
        parProp = par_i_old
        accRatio = 0
        muProp =  meanParOld_i
        varProp = varParOld_i
      }
      if(lfE != 0){
        parProp <- varProp <- muProp <- rep(NA, K)
        tmp <- arg/((1 +arg)^2)  # n K
        tmp1 <- arg/(1 +arg)  # n K
        varProp[freeEff] <- sapply(freeEff, function(x)
          ( ((alpha[x]^2)/(gt^2)) *sum( 1+sum( tmp[,x][ind[,x]])) )^(-1) ) # 1* length(freeEff)



        muProp[freeEff] <- sapply(freeEff, function(x)
          alpha[x]/gt *sum((Y_i[,x] -tmp1[,x])[ind[,x]])*varProp[x] +(par_i_old[x]) ) # 1*length(freeEff)

        parProp[freeEff] <- RcppTN::rtn(.mean = muProp[freeEff], .sd = sqrt(varProp[freeEff]), # 1*K
                                        .low = rep(-1, lfE), .high = rep(1, lfE), .checks = FALSE)

        parProp[fix_nullEff] <- par_i_old[fix_nullEff]
        muProp[fix_nullEff] <- meanParOld_i[fix_nullEff]
        varProp[fix_nullEff] <- varParOld_i[fix_nullEff]

        accRatio = sum(log(RcppTN::dtn(par_i_old[freeEff], .mean = muProp[freeEff], .sd = sqrt(varProp[freeEff]),
                                       .low = rep(-1, lfE), .high = rep(1, lfE),.checks = FALSE))) -
          sum(log(RcppTN::dtn(parProp[freeEff], .mean = meanParOld_i[freeEff], .sd = sqrt(varParOld_i[freeEff]),
                              .low = rep(-1, lfE), .high = rep(1, lfE), .checks = FALSE)))

      }
    }  # close var
  }
  return( list(parProp = parProp, muProp = muProp, varProp = varProp, accRatio = accRatio) )
}
