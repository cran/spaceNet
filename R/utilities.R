
simMultiNet <- function( n, K, D = 2, alphaRef = 1,
                         sender = NULL, receiver = NULL,
                         muA = 0, tauA = NULL, nuA = 3,
                         muB = 0, tauB = NULL, nuB = 3,
                         probabilities = FALSE, boundA = NULL){

  ##constants
  noSendRec <- is.null(sender) & is.null(receiver)
  if( is.null(boundA)){
    boundA <- if ( noSendRec ) log( (log(n))/(n -log(n)) ) else 0
  }
  if( !is.null(sender) & !is.null(receiver)) muA <- 0.5
  if ( is.null(tauA) ) tauA <- (K-1)/K
  if ( is.null(tauB) ) tauB <- (K-1)/K

  #hyperparameters
  varAlpha <- 1/rchisq( 1, df = nuA )
  muAlpha <- RcppTN::rtn( muA,sqrt( tauA*varAlpha ), boundA, Inf)

  varBeta <- 1/rchisq( 1, df = nuB )
  muBeta <- RcppTN::rtn( muB,  sqrt( tauB*varBeta ), 0, Inf)

  #standard parameters
  alpha <- RcppTN::rtn( rep( muAlpha, K), rep( sqrt( varAlpha), K ),
                        rep( boundA, K ), rep( Inf, K ) )
  alpha[1] <- alphaRef
  beta <- c( 1, RcppTN::rtn( rep( muBeta, K -1 ), rep( sqrt( varBeta), K -1 ),
                             rep( 0, K -1), rep( Inf, K -1 ) ))

  #sender & receiver
  theta <- gamma <- matrix( 0, ncol = K, nrow = n )
  if( !noSendRec){
    if( !is.null(sender ) ){
      if( sender == "const"){theta <- matrix( rep( runif( n, -1, 1 ), K ),
                                              ncol = K, byrow = FALSE )
      theta[1,] <- 1

      }
      if( sender == "var"){theta <- matrix( runif( n*K, -1, 1 ), ncol = K )
      theta[1,] <- 1}
    }
    if( !is.null(receiver ) ){
      if( receiver == "const"){gamma <- matrix( rep( runif( n, -1, 1 ), K ),
                                                ncol = K, byrow = FALSE )
      gamma[1,] <- 1}
      if( receiver == "var"){gamma <- matrix( runif( n*K, -1, 1 ), ncol = K )
      gamma[1,] <- 1}
    }
  }
  gt <- if ( is.null(sender) | is.null(receiver) ) 1 else 0.5
  if ( noSendRec ) gammaTheta <- 1 else gammaTheta <- list2array(
    lapply(1:K, function(k) matrix(gamma[,k], n,n, byrow = TRUE) + theta[,k]) )*gt

  #latent coordinates & distances
  z <- MASS::mvrnorm(n, rep(0, D ), diag(1, D))
  disz <- as.matrix( dist(z )^2 )

  probs <- list()
  for ( k in 1:K ){
    arg <- exp(alpha[k]*( if(noSendRec) 1 else gammaTheta[,,k]) -beta[k]*disz)
    diag(arg ) <- 0
    probs[[k]] <- arg/(1 +arg )
  }

  Y <- lapply(1:K, function(k)
    matrix( rbinom( n^2, 1, probs[[k]] ), n, n ) )

  parameters <- list( alpha = alpha, beta = beta, z = z,
                      gamma = if(is.null(receiver)) NULL else gamma,
                      theta = if(is.null(sender)) NULL else theta )
  hyperparameters <- list( varAlpha = varAlpha, muAlpha = muAlpha,
                           varBeta = varBeta, muBeta = muBeta)
  out <- list( Y = Y, parameters = parameters,
               hyperparameters = hyperparameters,
               probabilities = if( probabilities == FALSE) NULL else probs )
  return(out)
}



list2array <- function(L)
  # Convert a list of matrices to an array
{
  array( unlist(L), dim = c(dim(L[[1]]), length(L)) )
}
