#
#============== Functions to initialize multiNet
#


startZ <- function(Y, n, D, K)
  # Initialize latent positions and distances
{
  z <- vapply( 1:K, function(k) {
    distZ <- sna::geodist(Y[,,k])$gdist
    distZ[ is.infinite(distZ) ] <- 7^2
    cmdscale( as.dist(distZ), k = D )
  }, double(n*D) )
  z <- matrix(Rfast::rowmeans(z), ncol = D)

  z <- if ( D == 1 ) z/max(z) else apply(z, 2, function(u) u/max(u) )
  distZ <- 3 * Rfast::Dist(z, square = TRUE)

  # starting values for muZ and sigmaZ
  MeanZOld_i <- z
  VarZOld_i <- rep(0.1, n)

  out <- list( z = z, distZ = distZ, MeanZOld_i = MeanZOld_i, VarZOld_i = VarZOld_i )
  return(out)
}



startLogit <- function(Y, K, distZ, alpha1, noSendRec, ind)
  # Initialize logit parameters
{
  predZ <- -c(distZ)
  tmp <- vapply(1:K, function(k) {
    Rfast::glm_logistic(predZ[ind[,,k]], c( Y[,,k][ind[,,k]] ))$be[,1]
  }, double(2) )
  alpha <- tmp[1,] + 2
  if( any( is.na(alpha) ) ){ alpha <- rep(alpha1, K) }
  beta <- tmp[2,]
  if( any( is.na(beta) ) ){ beta <- rep(1, K) }
  alpha[1] <- alpha1
  alpha <- ifelse( abs(alpha) > 3, alpha / ( 10^nchar(abs(round(alpha))) ), alpha )
  if( noSendRec == FALSE ){ alpha[which(alpha < 0.001)] <- 0 }
  beta[2:K] <- beta[2:K] + abs( 1 - beta[1] )
  beta[1] <- 1
  beta <- ifelse(beta < 0, 0, beta)
  alpha[is.na(alpha)] <- 0.1
  beta[is.na(beta)] <-1

  # mean and variance paramters
  muAlpha <- mean(alpha)
  muBeta <- mean(beta)
  sigmaAlpha <- var(alpha) + 1
  sigmaBeta <- var(beta) + 1

  # get starting values for mu_alpha^k and sigma_alpha^k
  meanAlphaOld_k  <- alpha[2:K]
  varAlphaOld_k <- rep(1, K - 1)

  # get starting values for mu_beta^k and sigma_beta^k
  meanBetaOld_k <- beta[2:K]
  varBetaOld_k <- rep(1, K - 1)

  out <- list(alpha = alpha, muAlpha = muAlpha, sigmaAlpha = sigmaAlpha,
              meanAlphaOld_k = meanAlphaOld_k, varAlphaOld_k = varAlphaOld_k,
              beta = beta, muBeta = muBeta, sigmaBeta = sigmaBeta,
              meanBetaOld_k = meanBetaOld_k, varBetaOld_k = varBetaOld_k)
  return(out)
}



startSend <- function(Y, arcSumY, n, K, sender, theta0)
  # Initialize sender parameters
{
  if ( is.null(sender) ) {
    theta <- matrix(0, n,K)
    return(theta)
  } else {
    tmp <- Y
    tmp[which(tmp == 0)] <- 1
    tmp[is.na(tmp)] <- 0
    abs_i_k <- which( apply(tmp, 3, rowSums ) == 0 )

    theta <- switch(sender,
                    const = {
                      theta <- colSums( Y, na.rm = TRUE )
                      theta <- apply( theta,2, function(x) x/max(x) )
                      theta <- matrix(rowMeans(theta), ncol = K, nrow = n, byrow = FALSE)
                      theta <- (theta - 0.5)/0.5  # to range in -1 , 1
                      # MODIFICA
                      theta[which(theta != 0)] <- 0
                      theta[cbind(theta0, 1:K)] <- 1
                      theta <- ifelse( theta < -1, -1,  theta )
                      theta <- ifelse( theta > 1, 1,  theta )
                    },
                    var = {
                      theta <- colSums( aperm(Y, c(2,1,3)), na.rm = TRUE )
                      theta <- apply( theta,2, function(x) x/max(x))
                      theta <- (theta -0.5)/(0.5) # to range in -1, 1
                      # MODIFICA
                      theta[which(theta!=0)] <- 0
                      theta[cbind(theta0, 1:K)] <- 1
                      theta <- ifelse( theta < -1, -1, theta )
                      theta <- ifelse( theta > 1, 1, theta )
                    }
    )
    if( length(abs_i_k) > 0 & sender == "var" ){
      theta[abs_i_k ] <- NA
    }
  }

  return(theta)
}



startRec <- function(Y, arcSumY, n, K, receiver, gamma0)
  # Initialize receiver parameters
{
  if ( is.null(receiver) ) {
    gamma <- matrix(0, n,K)
    return(gamma)
  } else {
    tmp <- Y
    tmp[which( tmp == 0)] <- 1
    tmp[is.na(tmp)] <- 0
    abs_i_k = which(apply(tmp, 3, colSums ) == 0)
    gamma <- switch(receiver,
                    const = {
                      gamma <- colSums( Y, na.rm = TRUE )
                      gamma <- apply( gamma,2, function(x) x/max(x))
                      gamma <- matrix(rowMeans(gamma), ncol = K, nrow = n, byrow = F)
                      gamma <- (gamma -0.5)/0.5  # to range in -1 , 1
                      # MODIFICA
                      gamma[which(gamma!=0)] <- 0
                      gamma[cbind(gamma0, 1:K)] <- 1
                      gamma <- ifelse(gamma < -1, -1, gamma)
                      gamma <- ifelse(gamma > 1, 1, gamma)
                    },
                    var = {
                      gamma <- colSums( Y, na.rm = TRUE )
                      gamma <- apply( gamma,2, function(x) x/max(x))
                      gamma <- (gamma -0.5)/0.5  # to range in -1 , 1
                      # MODIFICA
                      gamma[which(gamma!=0)] <- 0
                      gamma[cbind(gamma0, 1:K)] <- 1
                      gamma <- ifelse(gamma < -1, -1, gamma)
                      gamma <- ifelse(gamma > 1, 1, gamma)
                    }
    )
    if( length(abs_i_k) > 0 & receiver == "var"){
      gamma[abs_i_k ] <- NA
    }
  }
  return(gamma)
}



startLambda <- function(covariates, nC)
  # Initialize parameter for covariates
{
  if ( !is.null(covariates) ) {
    # nC <- dim(covariates)[3]
    lambda <- rep(1, nC)
    meanLambda_l <- rep(0, nC)
    varLambda_l <- rep(1, nC)
    return( list(lambda = lambda, meanLambda_l = meanLambda_l, varLambda_l = varLambda_l) )
  } else {
    lambda <- 0
    return( list(lambda = lambda) )
  }
}



alphaRef <- function(Y, D = 2)
  # Intialize alpha in the reference network
{
  sender <- receiver <- NULL    # not implemented now
  if ( is.list(Y) ) Y <- list2array(Y)
  n <- nrow(Y[,, 1])
  pHat <- sum(Y[,,1], na.rm = TRUE)/(n*(n-1))
  out <- log( pHat/(1 -pHat) ) + 2
  if( !is.null(sender) | !is.null(receiver) & out <= 0 ) out <- 0.001
  return(out)
}
