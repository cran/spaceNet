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
                         noSendRec, const = FALSE)
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

  # term for acceptance ratio ---------------------
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

  if ( !const ) return( sum( llk + term ) )

  # return value with constant --------------------
  return( sum( llk + term ) +
            -0.5 *( K*( log(sigmaAlpha) +log(sigmaBeta) ) +
                      log(tauB*sigmaBeta) + log(tauA*sigmaAlpha) +
                      (muBeta - muB)^2/(tauB*sigmaBeta) + (muAlpha - muA)^2/(tauA*sigmaAlpha) +
                      (1/sigmaAlpha) + (1/sigmaBeta) ) -
            ( (nuB/2 + 1)*log(sigmaBeta) + (nuA/2 + 1)*log(sigmaAlpha) ) +
            if ( lambdaCovSum != 0 ) {
              -0.5*( sum( (muLambda_l - muL)^2/(tauL*sigmaLambda_l) ) + sum( log(tauL*sigmaLambda_l) ) +
                       sum( log(sigmaLambda_l) ) + sum(1/sigmaLambda_l) ) -
                sum( (nuL/2 + 1)*log(sigmaLambda_l) )
            } else 0
  )
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
