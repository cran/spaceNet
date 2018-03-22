#
#=========== MCMC for Latent Space Modeling of Multidimensional Networks
#

multiNet <- function( Y, niter = 1000, D = 2,
                      muA = 0, tauA = NULL, nuA = 3,
                      muB = 0, tauB = NULL, nuB = 3,
                      muL = 0, tauL = NULL, nuL = 3,
                      alphaRef = 0.1,
                      # sender = c("const", "var"), receiver = c("const", "var"),
                      covariates = NULL, DIC = FALSE,
                      burnIn = round(niter*0.3), trace = TRUE, allChains = FALSE,
                      refSpace = NULL )
{
  call <- match.call()

  # we like arrays
  if ( is.list(Y) ) Y <- list2array(Y)
  if ( !is.null(covariates) ) {
    if ( is.matrix(covariates) ){ covariates <- list(covariates)}
    if ( is.list(covariates) ) {
      covariates <- list2array(covariates)
    }
  }
  # no sender and receiver
  sender <- receiver <- NULL     # no sender and receiver for v1
  noSendRec <- is.null(sender) & is.null(receiver)

  # constants
  K <- dim(Y)[3]     # number of views
  n <- dim(Y)[1]     # number of nodes
  nn <- n*n
  # maximum correlation for 2 subsequential updates of the sender/receiver parameters
  boundGammaTheta <- 0.975
  boundA <- if ( noSendRec ) log( (log(n))/(n -log(n)) ) else 0   # lower bound for the alphas
  nC <- dim(covariates)[3]     # number of covariates, if any, else NULL
  if( !is.null(nC) ) {
    p <- if ( is.na(nC) ) 1 else seq(nC)
    if( is.na(nC) ) nC <- 1
  }

  # create indicator variable
  ind <- !is.na(Y)

  # store dBar values for DIC
  if( DIC ){
    dBar <- rep(NA, niter)
  }

  # tau hyperparameters
  if ( is.null(tauA) ) tauA <- (K-1)/K
  if ( is.null(tauB) ) tauB <- (K-1)/K
  if ( is.null(tauL) | !is.null(covariates) ){
    if ( is.null(nC) ) tauL <- 0.5
    if ( !is.null(nC) ){
      tauL <- (nC-1)/nC
      if ( is.na(tauL) | tauL == 0) tauL <- 0.5
    }
  }

  # sender/receiver
  if ( length(sender) == 2 ) sender <- NULL
  if ( length(receiver) == 2 ) receiver <- NULL
  if ( !is.null(sender) ) {
    sender <- match.arg(sender, c("const", "var"))
    theta0 <- Rfast::colMaxs( colSums(aperm(Y, c(2,1,3)), na.rm = TRUE) )
    if ( sender == "const" ){
      tmp <- Y
      tmp[which( tmp == 0)] <- 1
      tmp[is.na(tmp)] <- 0
      elegible <- rowSums( apply(tmp, 3, rowSums ) != 0) == K
      pick <- which( elegible[theta0] == TRUE )[1]
      theta0 <- rep(theta0[pick], K)
    }
  } else theta0 <- NULL
  if ( !is.null(receiver) ) {
    receiver <- match.arg(receiver, c("const", "var"))
    gamma0 <- Rfast::colMaxs( colSums(Y, na.rm = TRUE) )
    if( receiver == "const" ){
      tmp <- Y
      tmp[which( tmp == 0)] <- 1
      tmp[is.na(tmp)] <- 0
      elegible <- rowSums( apply(tmp, 3, colSums ) != 0) == K
      pick <- which( elegible[gamma0] == TRUE )[1]
      gamma0 <- rep(gamma0[pick], K)
    }
  } else gamma0 <- NULL
  noSendRec <- is.null(sender) & is.null(receiver)

  # number of arcs for each view
  arcSumY <- vapply( 1:K, function(k) sum(Y[,,k], na.rm = TRUE), double(1) )

  # starting values -------------------------------------------------------------------
  inZ <- startZ(Y, n, D, K)
  inLogit <- startLogit(Y, K, inZ$distZ, alphaRef, noSendRec, ind)
  theta <- startSend(Y, arcSumY, n, K, sender, theta0)
  gamma <- startRec(Y, arcSumY, n, K, receiver, gamma0)
  inLambda <- startLambda(covariates, nC)

  # objects to be stored --------------------------------------------------------------
  # nPar <- K - 1 + K - 1 + 4       # number of logit parameters
  ALPHA <- BETA <- accALPHABETA <- matrix(NA, niter + 1, K - 1)
  MSALPHA <- MSBETA <- accMSALPHA <- accMSBETA <-
    matrix(NA, niter + 1, 2)   # muAlpha and sigmaAlpha - muBeta and sigmaBeta
  ALPHA[1,] <- inLogit$alpha[2:K]
  BETA[1,] <- inLogit$beta[2:K]
  MSALPHA[1,] <- c(inLogit$muAlpha, inLogit$sigmaAlpha)
  MSBETA[1,] <- c(inLogit$muBeta, inLogit$sigmaBeta)
  Z <- array(NA, c(n, D, niter + 1))
  Z[,,1] <- inZ$z
  accZ <- matrix(NA, n, niter + 1)
  #
  meanThetaOld <- meanGammaOld <- varGammaOld <- varThetaOld <- NULL
  if ( !is.null(sender) ) {
    THETA <- accTHETA <- array(NA, c(n, K, niter + 1))
    THETA[,,1] <- theta
    meanThetaOld <- matrix(0, ncol = K, nrow = n)
    varThetaOld <- matrix(1, ncol = K, nrow = n)
  }
  if ( !is.null(receiver) ) {
    GAMMA <- accGAMMA <- array(NA, c(n, K, niter + 1))
    GAMMA[,,1] <- gamma
    meanGammaOld <- matrix(0, ncol = K, nrow = n)
    varGammaOld <- matrix( 1, ncol = K, nrow = n)
  }

  if ( !is.null(covariates) ) {
    # nC <- dim(covariates)[3]
    LAMBDA <- accLAMBDA  <- MLAMBDA <- SLAMBDA <- matrix(NA, niter + 1, nC)
  }

  # correlation with reference latent space -- for simulated data experiments
  corrZ <- if ( !is.null(refSpace) ) rep(NA, niter) else NULL


  # MCMC ------------------------------------------------------------------------------
  alpha <- inLogit$alpha[1:K]
  beta <- inLogit$beta[1:K]
  muAlpha <- inLogit$muAlpha                 # prior parameters
  muBeta <- inLogit$muBeta
  sigmaAlpha <- inLogit$sigmaAlpha
  sigmaBeta <- inLogit$sigmaBeta
  meanAlphaOld_k <- inLogit$meanAlphaOld_k   # proposal parameters
  varAlphaOld_k <- inLogit$varAlphaOld_k
  meanBetaOld_k <- inLogit$meanBetaOld_k
  varBetaOld_k <- inLogit$varBetaOld_k
  #
  distZ <- inZ$distZ
  z <- inZ$z
  meanZOld_i <- inZ$MeanZOld_i
  varZOld_i  <- inZ$VarZOld_i
  #
  if ( !is.null(covariates) ) {
    lambda <- inLambda$lambda
    muLambda_l <- inLambda$meanLambda_l     # prior parameters
    sigmaLambda_l <- inLambda$varLambda_l
    meanLambdaOld_l <- rep(0, nC)           # proposal parameters
    varLambdaOld_l <- rep(1, nC)
    #
    # initialize sufficient stats for lambda
    lambdaCov <- covariates * array( rep(lambda, each = nn), c(n, n, nC) )
    # faster than lambdaCovSum <- apply(lambdaCov, 1:2, sum)
    tmp1 <- c(lambdaCov)
    tmp2 <- paste0("tmp1[", (p-1)*(nn) + 1, ":", p*nn, "]", collapse = "+")
    lambdaCovSum <- matrix( eval(parse(text = tmp2)), n,n )
  } else {
    lambda <- inLambda$lambda
    lambdaCov <- NULL
    lambdaCovSum <- 0
  }

  gammaTheta <- 1

  pbar <- txtProgressBar(min = 2, max = (niter+1), style = 3)
  on.exit( close(pbar) )
  for ( it in 2:(niter + 1) ) {
    # if ( it%%10 == 0 ) showTrace(it)
    setTxtProgressBar(pbar, it)

    # sufficient stats for gamma and theta (only if they are in the model)
    gt <- if ( is.null(sender) | is.null(receiver) ) 1 else 0.5
    if ( !noSendRec ) {
      gammaTheta <- list2array(
        lapply(1:K, function(k) matrix(gamma[,k], n,n, byrow = TRUE) + theta[,k]) )*gt
    }

    ### alpha and beta parameters ...................................
    # update sigmaAlpha
    sigmaAlpha <- sigmaFc(alpha, muAlpha, nuA, tauA) * (n/5)
    MSALPHA[it,2] <- sigmaAlpha
    # update muAlpha
    muAlpha <- muFc(alpha, sigmaAlpha, tauA, muA, boundA, Inf)
    MSALPHA[it,1] <- muAlpha

    # update sigmaBeta
    sigmaBeta <- sigmaFc(beta, muBeta, nuB, tauB)
    MSBETA[it,2] <- sigmaBeta
    # update muBeta
    muBeta <- muFc(beta, sigmaBeta, tauB, muB, 0, Inf)
    MSBETA[it,1] <- muBeta

    if ( !is.null(covariates) ) {
      # update sigmaLambda_l
      sigmaLambda_l <- vapply( 1:nC, function(l) sigmaFc(lambda[l], muLambda_l[l], nuL, tauL), numeric(1) )
      SLAMBDA[it,] <- sigmaLambda_l
      # update muLambda_l
      muLambda_l <- vapply( 1:nC, function(l) muFc(lambda[l], sigmaLambda_l[l], tauL, muL, 0, Inf), numeric(1) )
      MLAMBDA[it,] <- muLambda_l
    }

    # update alpha
    alphaBetaCand <-
      vapply( 2:K, function(k) {
        alphaBetaProp_k(k, Y[,,k], n, arcSumY[k], alpha[k], beta[k],
                        if (noSendRec) gammaTheta else gammaTheta[,,k], distZ,
                        muAlpha, sigmaAlpha, muBeta, sigmaBeta,
                        boundA, lambdaCovSum, ind[,,k], noSendRec)
      }, numeric(6) )
    alphaTry <- alphaBetaCand[1,]
    alphaTryMean <- alphaBetaCand[2,]
    alphaTryVar <- alphaBetaCand[3,]
    betaTry <- alphaBetaCand[4,]
    betaTryMean <- alphaBetaCand[5,]
    betaTryVar <- alphaBetaCand[6,]
    #
    alphaOld <- alpha     # keep previous values
    alphaNew <- alphaOld
    betaOld <- beta
    betaNew <- betaOld

    # compute acceptance ratio
    for ( k in 2:K ) {
      # out <- vapply(2:K, function(k) {
      alphaNew[k] <- alphaTry[k-1]
      betaNew[k] <- betaTry[k-1]

      logDiff <- logPosterior(Y, n, K, alphaNew, betaNew,
                              muAlpha, sigmaAlpha, muBeta, sigmaBeta,
                              muLambda_l, sigmaLambda_l,
                              muA, tauA, nuA,
                              muB, tauB, nuB,
                              muL, tauL, nuL,
                              gammaTheta, 0, lambda,
                              z, distZ, lambdaCovSum,
                              ind, term = "alphaBeta", noSendRec) -
        logPosterior(Y, n, K, alphaOld, betaOld,
                     muAlpha, sigmaAlpha, muBeta, sigmaBeta,
                     muLambda_l, sigmaLambda_l,
                     muA, tauA, nuA,
                     muB, tauB, nuB,
                     muL, tauL, nuL,
                     gammaTheta, 0, lambda,
                     z, distZ, lambdaCovSum,
                     ind, term = "alphaBeta", noSendRec)

      num <- c( RcppTN::dtn(alphaOld[k], .mean = alphaTryMean[k-1], .sd = sqrt(alphaTryVar[k-1]),
                            .low = boundA, .checks = FALSE),
                RcppTN::dtn(betaOld[k], .mean = betaTryMean[k-1], .sd = sqrt(betaTryVar[k-1]),
                            .low = 0, .checks = FALSE)
      )
      den <- c( RcppTN::dtn(alphaNew[k], .mean = meanAlphaOld_k[k-1], .sd = sqrt(varAlphaOld_k[k-1]),
                            .low = boundA, .checks = FALSE),
                RcppTN::dtn(betaNew[k], .mean = meanBetaOld_k[k-1], .sd = sqrt(varBetaOld_k[k-1]),
                            .low = 0, .checks = FALSE)
      )

      logAccAlphaBeta <- logDiff + ( sum( log(num) ) - sum( log(den) ) )
      logU <- log( runif(1) )
      if ( logAccAlphaBeta < logU | any(is.nan(num)) | any(is.infinite(num)) | any(num == 0) ) {
        # reject
        alphaNew[k] <- alphaOld[k]
        betaNew[k] <- betaOld[k]
        # accAlphaBeta <- 0
        accALPHABETA[it,k-1] <- 0
      } else {
        # accept
        alphaOld[k] <- alphaNew[k]
        meanAlphaOld_k[k-1] <- alphaTryMean[k-1]
        varAlphaOld_k[k-1] <- alphaTryVar[k-1]
        betaOld[k] <- betaNew[k]
        meanBetaOld_k[k-1] <- betaTryMean[k-1]
        varBetaOld_k[k-1] <- betaTryVar[k-1]
        # accAlphaBeta <- 1
        accALPHABETA[it,k-1] <- 1
      }
      # return( matrix( c(alphaNew[k], meanAlphaOld_k[k-1], varAlphaOld_k[k-1],
      #                   betaNew[k], meanBetaOld_k[k-1], varBetaOld_k[k-1], accAlphaBeta), 7, 1 ) )
      # }, numeric(7))
    }

    alpha <- alphaOld
    beta <- betaOld
    #
    ALPHA[it,] <- alpha[-1]
    BETA[it,] <- beta[-1]
    ###..............................................................

    ###..............................................................
    zOld <- z
    distZOld <- distZ
    for ( i in 1:n ) {
      zProp <- zProp_i(i, Y[i,,], n, K, D, z, distZ[i,], alpha, beta,
                       ind[i,,],  if (noSendRec) gammaTheta else gammaTheta[i,,],
                       if ( length(lambdaCovSum) != 1 ) lambdaCovSum[i,] else 0)

      zNew <- z
      zNew[i,] <- zProp$zNew_i
      distZNew <- Rfast::Dist(zNew, square = TRUE)

      logDiff <- logPosterior(Y, n, K, alpha, beta,
                              muAlpha, sigmaAlpha,
                              muBeta, sigmaBeta,
                              muLambda_l, sigmaLambda_l,
                              muA, tauA, nuA,
                              muB, tauB, nuB,
                              muL, tauL, nuL,
                              gammaTheta, 0, lambda,
                              zNew, distZNew,
                              lambdaCovSum,
                              ind, term = "z", noSendRec) -
        logPosterior(Y, n, K, alpha, beta,
                     muAlpha, sigmaAlpha,
                     muBeta, sigmaBeta,
                     muLambda_l, sigmaLambda_l,
                     muA, tauA, nuA,
                     muB, tauB, nuB,
                     muL, tauL, nuL,
                     gammaTheta, 0, lambda,
                     z, distZ,
                     lambdaCovSum,
                     ind, term = "z", noSendRec)

      logAccZ <- logDiff +
        ( Rfast::dmvnorm( z[i,], zProp$meanZ_i, diag(zProp$varZ_i, D), logged = TRUE ) -
            Rfast::dmvnorm( zNew[i,], meanZOld_i[i,], diag(varZOld_i[i], D), logged = TRUE )
        )
      logU <- log( runif(1) )

      if ( logAccZ >= logU ) {
        z[i,] <- zProp$zNew_i
        distZ <- distZNew
        meanZOld_i[i,] <- zProp$meanZ_i
        varZOld_i[i] <- zProp$varZ_i
        accZ[i,it] <- 1
      } else accZ[i,it] <- 0
    }

    ##--- check if the new set is just a rotation of the previous-----
    rotCheck <- vegan::protest( z, zOld, scale = FALSE, translation = FALSE,
                                permutations = permute::how(nperm = 99) )[[6]]
    if ( rotCheck >= 0.95 ) {
      z <- zOld
      distZ <- distZOld
    }
    Z[,,it] <- z
    ##-----------------------------------------------------------------

    ##------correlation with a ref latent space------------------------
    if ( !is.null(refSpace) ){
      corrZ[it] <- vegan::protest( z, refSpace, scale = FALSE, translation = FALSE )[[6]]
    }
    ##-----------------------------------------------------------------
    #close update only z

    ###..............................................................
    if( !noSendRec){
      if(!is.null(sender)){
        thetaOld <- theta
        gammaThetaOld <- gammaTheta
        for ( i in 1:n ) {
          ##-----propose a new value for theta_i ---------------------------------
          thetaProp <- gammaThetaProp_i(theta[i,], Y[i,,], n, K, distZ[i,],i,
                                        alpha, beta, gt, gammaTheta[i,,], gamma,
                                        if ( length(lambdaCovSum) != 1 ) lambdaCovSum[i,] else 0,
                                        sender, theta0,
                                        if ( !is.null(sender) ) meanThetaOld[i,] else 0,
                                        if ( !is.null(sender) ) varThetaOld[i,] else 0,
                                        ind[i,,]
          )  ##effect is either N or C or V
          thetaNew <- theta
          thetaNew[i, ] <- thetaProp$parProp

          gammaThetaNew <- list2array(
            lapply(1:K, function(k) matrix(gamma[,k], n,n, byrow = TRUE) + thetaNew[,k]) )*gt

          logDiff <- logPosterior(Y, n, K, alpha, beta,
                                  muAlpha, sigmaAlpha,
                                  muBeta, sigmaBeta,
                                  muLambda_l, sigmaLambda_l,
                                  muA, tauA, nuA,
                                  muB, tauB, nuB,
                                  muL, tauL, nuL,
                                  gammaThetaNew,0, lambda,
                                  z, distZ,
                                  lambdaCovSum,
                                  ind, term = "sendRec", noSendRec) -
            logPosterior(Y, n, K, alpha, beta,
                         muAlpha, sigmaAlpha,
                         muBeta, sigmaBeta,
                         muLambda_l, sigmaLambda_l,
                         muA, tauA, nuA,
                         muB, tauB, nuB,
                         muL, tauL, nuL,
                         gammaTheta, 0, lambda,
                         z, distZ,
                         lambdaCovSum,
                         ind, term ="sendRec", noSendRec)

          logAccTheta <-  logDiff + thetaProp$accRatio

          logU <- log( runif(1) )

          if ( logAccTheta >= logU & !is.na(logAccTheta) ) {
            gammaTheta <- gammaThetaNew
            theta[i, ] <- thetaNew[i, ]
            meanThetaOld[i,] <- thetaProp$muProp
            varThetaOld[i,] <- thetaProp$varProp
            accTHETA[i,,it] <- rep(1, K)
          } else{
            accTHETA[i,,it] <- rep(0, K)
          }
        } # end i loop

        if ( !is.null(receiver) ) {
          checkTheta <- sapply( 1:K, function(x)
            cor( c(alphaRef, ALPHA[it-1,])[x]*c(gammaThetaOld[,,x]),
                 alpha[x]*c(gammaTheta[,,x]), use = "na.or.complete" )
          ) }
        else{ checkTheta <- sapply(1:K, function(x)
          cor( c(alphaRef, ALPHA[it-1,])[x] *c(thetaOld[,x]),
               alpha[x]*c(theta[,x]), use = "na.or.complete" )
        ) }
        subs <- which( abs(checkTheta) - boundGammaTheta > 0)

        if(length( subs) >0){
          theta[,subs] <- thetaOld[,subs]
          gammaTheta[,,subs] <- gammaThetaOld[,,subs]
        }
        THETA[,,it] <- theta
      }


      if(!is.null(receiver)){
        gammaOld <- gamma
        gammaThetaOld <- gammaTheta
        for ( i in 1:n ) {
          ##-----propose a new value for gamma_i ---------------------------------
          gammaProp <- gammaThetaProp_i(gamma[i,], Y[,i,], n, K, distZ[i,],i, # gamma
                                        alpha, beta, gt, gammaTheta[,i,], theta,
                                        if ( length(lambdaCovSum) != 1 ) lambdaCovSum[,i] else 0,
                                        receiver, gamma0,
                                        if ( !is.null(receiver) ) meanGammaOld[i,] else 0,
                                        if ( !is.null(receiver) ) varGammaOld[i,] else 0,
                                        ind[,i,]
          )  # effect is either N or C or V
          gammaNew <- gamma
          gammaNew[i, ] <- gammaProp$parProp

          gammaThetaNew <- list2array(
            lapply(1:K, function(k) matrix(gammaNew[,k], n,n, byrow = TRUE) + theta[,k]) )*gt

          logDiff <- logPosterior(Y, n, K, alpha, beta,
                                  muAlpha, sigmaAlpha,
                                  muBeta, sigmaBeta,
                                  muLambda_l, sigmaLambda_l,
                                  muA, tauA, nuA,
                                  muB, tauB, nuB,
                                  muL, tauL, nuL,
                                  gammaThetaNew,0, lambda,
                                  z, distZ,
                                  lambdaCovSum,
                                  ind, term = "sendRec", noSendRec) -
            logPosterior(Y, n, K, alpha, beta,
                         muAlpha, sigmaAlpha,
                         muBeta, sigmaBeta,
                         muLambda_l, sigmaLambda_l,
                         muA, tauA, nuA,
                         muB, tauB, nuB,
                         muL, tauL, nuL,
                         gammaTheta, 0, lambda,
                         z, distZ,
                         lambdaCovSum,
                         ind, term ="sendRec", noSendRec)

          logAccGamma <-  logDiff + gammaProp$accRatio
          logU <- log( runif(1) )

          if ( logAccGamma >= logU & !is.na(logAccGamma) ) {
            gammaTheta <- gammaThetaNew
            gamma[i, ] <- gammaNew[i, ]
            meanGammaOld[i,] <- gammaProp$muProp
            varGammaOld[i,] <- gammaProp$varProp
            accGAMMA[i,,it] <- rep(1, K)
          } else{
            accGAMMA[i,,it] <- rep(0, K)
          }
        } # end i loop

        if ( !is.null(sender) ) {
          checkGamma <- sapply(1:K, function(x)
            cor( c(alphaRef, ALPHA[it-1,])[x]*c(gammaThetaOld[,,x]),
                 alpha[x]*c(gammaTheta[,,x]), use = "na.or.complete" ))
        } else { checkGamma <- sapply(1:K, function(x)
          cor( c(alphaRef, ALPHA[it-1,])[x] *c(gammaOld[,x]),
               alpha[x]*c(gamma[,x]), use = "na.or.complete" ))
        }
        subs <- which( abs(checkGamma) - boundGammaTheta > 0)

        if( length(subs) > 0 ){
          gamma[,subs] <- gammaOld[,subs]
          gammaTheta[,,subs] <- gammaThetaOld[,,subs]
        }
        GAMMA[,,it] <- gamma
      }
    } # close update sender and receiver
    ###..............................................................

    ### lambda ......................................................
    if ( !is.null(covariates) ) {
      out <- vapply(1:nC, function(l) {
        lambdaProp_l(l, Y, n, K, nC, lambdaCov, lambdaCovSum, muLambda_l, sigmaLambda_l,
                     alpha, beta, noSendRec, gammaTheta, distZ, ind)
      }, numeric(3))

      lambdaTry <- out[1,]
      lambdaTryMean <- out[2,]
      lambdaTryVar <- out[3,]
      lambdaOld <- lambda
      lambdaNew <- lambdaOld
      lambdaCovOld <- lambdaCov
      lambdaCovSumOld <- lambdaCovSum


      for( l in 1:nC ){
        lambdaNew[l] <- lambdaTry[l]

        # update sufficient statistics
        lambdaCovNew <- covariates * array( rep(lambdaNew, each = nn), c(n, n, nC) )
        tmp1 <- c(lambdaCovNew)
        tmp2 <- paste0("tmp1[", (p-1)*(nn) + 1, ":", p*nn, "]", collapse = "+")
        lambdaCovSumNew <- matrix( eval(parse(text = tmp2)), n,n )

        logDiff <- logPosterior(Y, n, K, alpha, beta,
                                muAlpha, sigmaAlpha, muBeta, sigmaBeta,
                                muLambda_l, sigmaLambda_l,
                                muA, tauA, nuA,
                                muB, tauB, nuB,
                                muL, tauL, nuL,
                                gammaTheta, 0, lambdaNew,
                                z, distZ, lambdaCovSumNew,
                                ind, term = "lambda", noSendRec) -
          logPosterior(Y, n, K, alpha, beta,
                       muAlpha, sigmaAlpha, muBeta, sigmaBeta,
                       muLambda_l, sigmaLambda_l,
                       muA, tauA, nuA,
                       muB, tauB, nuB,
                       muL, tauL, nuL,
                       gammaTheta, 0, lambdaOld,
                       z, distZ, lambdaCovSumOld,
                       ind, term = "lambda", noSendRec)

        num <- RcppTN::dtn(lambdaOld[l], .mean = lambdaTryMean[l], .sd = sqrt(lambdaTryVar[l]),
                           .low = 0, .checks = FALSE)
        den <- RcppTN::dtn(lambdaNew[l], .mean = meanLambdaOld_l[l], .sd = sqrt(varLambdaOld_l[l]),
                           .low = 0, .checks = FALSE)

        logLambda <- logDiff + ( log(num) - log(den) )
        logU <- log( runif(1) )
        if ( logLambda < logU | any(is.nan(num)) | any(is.infinite(num)) | any(num == 0) ) {
          # reject
          lambdaNew[l] <- lambdaOld[l]
          accLAMBDA[it,l] <- 0
        } else {
          # accept
          lambdaOld[l] <- lambdaNew[l]
          meanLambdaOld_l[l] <- lambdaTryMean[l]
          varLambdaOld_l[l] <- lambdaTryVar[l]
          accLAMBDA[it,l] <- 1
          lambdaCovSumOld <- lambdaCovSumNew
          lambdaCovOld <- lambdaCovNew
        }
      }
      lambdaCovSum <- lambdaCovSumOld
      lambda <- lambdaOld
      lambdaCov <- lambdaCovOld
      LAMBDA[it,] <- lambda
    }
    ###..............................................................

    # compute first part DIC
    if( DIC ){
      dBar[it-1] <- -2*loglik(Y, n, K, alpha, beta,
                              gammaTheta, distZ,
                              lambdaCovSum,ind,
                              noSendRec)
    }
  } # mcmc loop

  # output of the chain
  set <- (burnIn + 2):(niter + 1)
  ALPHA_mean <- colMeans( ALPHA[set, ])
  ALPHA_sd <- apply( ALPHA[set, ], 2, sd)
  BETA_mean <- colMeans( BETA[set, ])
  BETA_sd <- apply( BETA[set, ], 2, sd)
  accALPHABETA_rate <- colSums(accALPHABETA[set,])/(niter -burnIn -2)
  #### accALPHABETA

  accZ_rate <- rowSums( accZ[,set] )/(niter -burnIn -2)
  Z_mean <- apply(Z[,,set], c(1,2), mean)
  Z_sd <- apply(Z[,,set], c(1,2), sd)
  distZ_mean <- as.matrix( dist( Z_mean )^2 )

  if( !is.null(receiver) ) {
    GAMMA_mean <- apply( GAMMA[,,set], c(1,2), mean)
    GAMMA_sd <- apply( GAMMA[,,set], c(1,2), sd)
    accGAMMA_rate <- apply(accGAMMA[,,set], c(1,2), sum)/(niter -burnIn -2)
  } else GAMMA_mean <- gamma
  if( !is.null(sender) ){
    THETA_mean <- apply( THETA[,,set], c(1,2), mean)
    THETA_sd <- apply( THETA[,,set], c(1,2), sd)
    accTHETA_rate <- apply( accTHETA[,,set], c(1,2), sum )/(niter -burnIn -2)
  } else THETA_mean <- theta
  if( !noSendRec ) {
    gammaTheta_mean <- list2array(
      lapply(1:K, function(k) matrix(GAMMA_mean[,k], n,n, byrow = TRUE) + THETA_mean[,k]) )*gt
  } else gammaTheta_mean <- 1

  if( !is.null(covariates) ){
    accLAMBDA_rate <- colSums(accLAMBDA[set,,drop = FALSE])/(niter -burnIn -2)
    LAMBDA_mean <- colMeans(LAMBDA[set,,drop = FALSE])
    LAMBDA_sd <- apply( ALPHA[set,,drop = FALSE], 2, sd )
    lambdaCov_mean <- covariates * array( rep(LAMBDA_mean, each = nn), c(n, n, nC) )
    tmp1 <- c(lambdaCov_mean)
    tmp2 <- paste0("tmp1[", (p-1)*(nn) + 1, ":", p*nn, "]", collapse = "+")
    lambdaCovSum_mean <- matrix( eval(parse(text = tmp2)), n,n )
  } else lambdaCovSum_mean <- lambdaCovSum


  # todo: lambda cov sum mean
  if ( DIC ){
    dHat <- -2*loglik(Y, n, K, c(alphaRef, ALPHA_mean), c(1, BETA_mean),
                      gammaTheta_mean, distZ_mean,
                      lambdaCovSum_mean, ind,
                      noSendRec)
    dBar <- mean( dBar[(burnIn+1):niter] )
    DIC <- dHat + 2*( dBar - dHat )
  }

  # prepare output
  parameters <- list(alpha = list(mean = c(alphaRef, ALPHA_mean), sd = c(0, ALPHA_sd)),
                     beta = list(mean = c(1, BETA_mean), sd = c(0, BETA_sd)),
                     gamma = if ( !is.null(receiver) ) {
                       list(mean = GAMMA_mean, sd = GAMMA_sd) } else (list(mean = NULL, sd = NULL)),
                     theta = if ( !is.null(sender) ) {
                       list(mean = THETA_mean, sd = THETA_sd) } else (list(mean = NULL, sd = NULL)),
                     lambda = if ( !is.null(covariates) ) {
                       list(mean = LAMBDA_mean, sd = LAMBDA_sd) } else (list(mean = NULL, sd = NULL))
  )
  latPos <- list(mean = Z_mean, sd = Z_sd)
  accRates <- list(alpha = accALPHABETA_rate[1], beta = accALPHABETA_rate[2],
                   gamma = if (!is.null(receiver)) accGAMMA_rate else NULL,
                   theta = if (!is.null(sender)) accTHETA_rate else NULL,
                   lambda = if (!is.null(covariates)) accLAMBDA_rate else NULL,
                   latPos = accZ_rate)
  if ( allChains ) {
    allChains <- list(parameters = list(
      alpha = cbind(alphaRef, ALPHA[2:(niter+1),]), beta = cbind(1, BETA[2:(niter+1),]),
      gamma = if (!is.null(receiver)) GAMMA[,,2:(niter+1)] else NULL,
      theta = if (!is.null(sender)) THETA[,,2:(niter+1)] else NULL,
      lambda = if (!is.null(covariates)) LAMBDA[2:(niter+1),] else NULL),
      latPos = Z[,,2:(niter+1)],
      priorParameters = list(
        alphaPrior = MSALPHA, betaPrior = MSBETA,
        lambdaPrior = if (!is.null(covariates)) list(meanLambda = MLAMBDA, sdLambda = SLAMBDA) else NULL
      )
      # , accept = list(alpha = accALPHABETA[2:(niter+1),1], beta = accALPHABETA[2:(niter+1),2],
      #               gamma = if (!is.null(receiver)) accGAMMA[,,2:(niter+1)] else NULL,
      #               theta = if (!is.null(sender)) THETA[,,2:(niter+1)] else NULL,
      #               lambda = if (!is.null(covariates)) LAMBDA[2:(niter+1),] else NULL,
      #               latPos = accZ[,2:(niter+1)])
    )
  } else allChains <- NULL

  if ( !is.numeric(DIC) ) DIC <- NULL
  out <- list(n = n, K = K, D = D, parameters = parameters,
              latPos = latPos, accRates = accRates,
              DIC = DIC, allChains = allChains, corrRefSpace = corrZ,
              info = list(call = call, niter = niter, burnIn = burnIn,
                          # receiver = receiver, sender = sender,
                          covariates = if( !is.null(covariates) ) TRUE else FALSE,
                          L = nC)
  )
  class(out) <- "multiNet"
  return(out)
}


print.multiNet <- function(x, ...)
{
  alpha <- x$parameters$alpha$mean
  beta <- x$parameters$beta$mean
  # gamma <- x$parameters$gamma$mean
  # theta <- x$parameters$theta$mean
  lambda <- x$parameters$lambda$mean
  dic <- x$DIC
  nC <- x$info$L

  if ( is.null(nC) ) {
    h1 <- paste(x$D, "dimensional latent space model")
    h2 <- "    for multivariate networks"
    sep <- paste0( rep("=", max(nchar(h1), nchar(h2)) + 5), collapse = "" )
    cat("\n", sep, "\n")
    cat("  ", h1, "\n")
    cat("  ", h2, "\n")
    cat("", sep, "\n", "\n")
    cat(" View 1", " --- ", "P(y = 1) =", round(alpha[1], 2), "- 1.00*dZ", "\n")
    for ( k in 2:x$K ) {
      cat(" View", k, " --- ", "P(y = 1) =", round(alpha[k], 2), "-",
          paste(round(beta[k], 2), "dZ", sep = "*"), "\n")
    }
  } else {
    covs <- paste0("X", 1:nC)
    h1 <- paste("    ", x$D, "dimensional latent space model")
    h2 <- "for multivariate networks, with covariates"
    sep <- paste0( rep("=", max(nchar(h1), nchar(h2)) + 5), collapse = "" )
    cat("\n", sep, "\n")
    cat("  ", h1, "\n")
    cat("  ", h2, "\n")
    cat("", sep, "\n", "\n")
    cat(" View 1", " --- ", "P(y = 1) =", round(alpha[1], 2), "- 1.00*dZ", "-",
        paste(paste(round(lambda, 2), covs, sep = "*"), collapse = " - "), "\n" )
    for ( k in 2:x$K ) {
      cat(" View", k, " --- ", "P(y = 1) =", round(alpha[k], 2), "-",
          paste(round(beta[k], 2), "dZ", sep = "*"),"-",
          paste(paste(round(lambda, 2), covs, sep = "*"), collapse = " - "),
          "\n")
    }
  }
  if ( !is.null(dic) ) cat("\n", "DIC = ", dic, "\n")
  cat("\n")
}
