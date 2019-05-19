#
#============== Functions for computing log-posterior distributions and loglik
#


logPosterior <- function(Y, n, K, alpha, beta,
                         muAlpha, sigmaAlpha,
                         muBeta, sigmaBeta,
                         muLambda_l, sigmaLambda_l,
                         muA, tauA, nuA,
                         muB, tauB, nuB,
                         muL, tauL, nuL,
                         gammaTheta, effect,
                         lambda,
                         z, distZ,
                         lambdaCovSum, ind,
                         term = c("alphaBeta", "z", "lambda", "sendRec"),
                         noSendRec, individual = FALSE)
# Compute model log-posterior
{
  # likelihood term -------------------------------
  if( individual == FALSE){
    llk <-
      vapply(1:K, function(k) {
        arg <- alpha[k]* ( if(noSendRec) gammaTheta else gammaTheta[,,k] ) - beta[k]*distZ - lambdaCovSum
        tmp <- log( 1 + exp(arg) )
        diag(tmp) <- 0
        sum( (Y[,,k]*arg - tmp)[ind[,,k]] )
      }, numeric(1) )
  } else {
    llk <-
      vapply(1:K, function(k) {
        arg <- alpha[k]* ( if(noSendRec) gammaTheta else gammaTheta[,,k] ) - beta[k]*distZ - lambdaCovSum
        tmp <- log( 1 + exp(arg) )
        diag(tmp) <- 0
        (Y[,,k]*arg - tmp ) ##diagonal is included
      }, numeric(n^2) )
    llk <- rowSums( llk)

    alphaBeta1 = sum( vapply(1:K, function(k) {
      -0.5*( (alpha[k] - muAlpha)^2/sigmaAlpha + (beta[k] - muBeta)^2/sigmaBeta )},
      numeric(1) ))

    zTemp <- -0.5*rowSums(z^2)
    zPr <- c( zTemp %*%t(zTemp) )

    lambdap <- if( lambdaCovSum != 0 ) -0.5*sum( (lambda - muLambda_l)^2/sigmaLambda_l ) else 0

    outTemp <- llk + alphaBeta1 + zPr +lambdap
    out <- outTemp[ -seq(1, n^2, by = n) ]
  }
  # term for acceptance ratio ---------------------
  if( individual == FALSE){
    term <- switch(term,
                   alphaBeta = {
                     vapply(1:K, function(k) {
                       -0.5*( (alpha[k] - muAlpha)^2/sigmaAlpha + (beta[k] - muBeta)^2/sigmaBeta )},
                       numeric(1) )
                   },
                   z = {
                     -0.5*sum(z^2)
                   },
                   lambda = {
                     -0.5*sum( (lambda - muLambda_l)^2/sigmaLambda_l )
                   },
                   sendRec = {
                     0   # uniform prior on the sender and receiver parameters
                   })
    out <- sum( llk + term )
  }
  # return value with constant --------------------
  return( out )
}



loglik <- function(Y, n, K, alpha, beta,
                   gammaTheta, distZ,
                   lambdaCovSum, ind,
                   noSendRec)
  # Compute model log-posterior
{
  # likelihood term -------------------------------
  llk <-
    vapply(1:K, function(k) {
      arg <- alpha[k]* ( if(noSendRec) gammaTheta else gammaTheta[,,k] ) - beta[k]*distZ - lambdaCovSum
      tmp <- log( 1 + exp(arg) )
      diag(tmp) <- 0
      sum( (Y[,,k]*arg - tmp)[ind[,,k]] )
    }, numeric(1) )

  return(sum(llk))
}
