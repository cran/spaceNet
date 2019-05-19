heuSearch <- function( Y, D = 2, nrep = 1000){
  # observed data
  DATA <- list2array(Y)
  n <- dim(DATA)[1]; K <- dim(DATA)[3]
  alphaRef <- alphaRef(DATA, D)

  # constants
  nmod = 9
  # storing
  temp <- matrix( NA, nrow = nrep*nmod, ncol = 7) # store the indexes
  colnames(temp) <- c("mCs_sdCs", "mCr_sdCr", "mCsT", "sdCsT",
                      "mCrT", "sdCrT", "mod")
  # computing the indexes
  index = 1
  cat("\n Simulating networks \n")
  for( t in 1:nrep){
    for ( m in 1:nmod){
      if( m == 1 ){ sender <- receiver <- NULL}
      if( m == 2 ){ receiver <- NULL ; sender <- "const" }
      if( m == 3 ){ sender <- NULL ; receiver <- "const" }
      if( m == 4 ){ sender <- receiver <- "const"}
      if( m == 5 ){ receiver <- NULL ; sender <- "var" }
      if( m == 6 ){ sender <- NULL ; receiver <- "var" }
      if( m == 7 ){ receiver <- "const" ; sender <- "var" }
      if( m == 8 ){ receiver <- "var" ; sender <- "const" }
      if( m == 9 ){ sender <- receiver <- "var" }

      data <- simMultiNet(n, K, sender = sender,
                          receiver = receiver, boundA = 0,
                          alphaRef = alphaRef, muA = 0.5)
      Yt <- list2array(data$Y)

      indegs <- sapply(1:K, function(x) colSums(Yt[,,x]))
      indegs <- indegs/max(indegs)
      outdegs <- sapply(1:K, function(x) rowSums(Yt[,,x]))
      outdegs <- outdegs/max(outdegs)
      c1 <- suppressWarnings( cor(indegs, use = "na.or.complete") ) #rk
      c2 <- suppressWarnings( cor(outdegs, use = "na.or.complete") ) #sk
      c3 <- suppressWarnings( cor(t(indegs), use = "na.or.complete") ) #ri
      c4 <- suppressWarnings( cor(t(outdegs), use = "na.or.complete") ) #si

      t1 <- mean(c(c2),  na.rm=TRUE) #sender K
      t2 <- sd(c(c2),  na.rm=TRUE) #sender K
      temp[index,1] <- abs( t1 -t2)
      t1 <-  mean(c(c1), na.rm=TRUE) #receiver K
      t2 <- sd(c(c1), na.rm=TRUE) #receiver K
      temp[index,2] <-  abs( t1 -t2)
      temp[index,3] <- mean(c(c4),  na.rm=TRUE) #sender i
      temp[index,4] <- sd(c(c4),  na.rm=TRUE) #sender i
      temp[index,5] <- mean(c(c3),  na.rm=TRUE) #receiver i
      temp[index,6] <- sd(c(c3),  na.rm=TRUE) #receiver i
      temp[index,7] <- m
      index = index +1
    } #close mod
  } # close nrep

  # DA
  cat("\n Heuristic search \n", "\n")
 browser()
  mclass <- factor(temp[,7])
  mp <- mclust::MclustDA( temp[, -7], class = mclass, verbose = FALSE)
  pred_mp <- mclust::predict.MclustDA(mp)
  out_mod_class <- round( table(predicted = pred_mp$classification, class = temp[,7])/ nrep, 3 )
  namesMod <- c("NN", "CN", "NC", "CC","VN", "NV", "VC","CV", "VV")
  colnames(out_mod_class) <- namesMod
  rownames(out_mod_class) <- namesMod[apply(out_mod_class,2, which.max)]

  ## compute class-model probabilities for observed data
  indegs <- sapply(1:K, function(x) colSums(DATA[,,x]))
  indegs <- indegs/max(indegs)
  outdegs <- sapply(1:K, function(x) rowSums(DATA[,,x]))
  outdegs <- outdegs/max(outdegs)
  c1 <- suppressWarnings( cor(indegs, use = "na.or.complete") ) #rk
  c2 <- suppressWarnings( cor(outdegs, use = "na.or.complete") ) #sk
  c3 <- suppressWarnings( cor(t(indegs), use = "na.or.complete") ) #ri
  c4 <- suppressWarnings( cor(t(outdegs), use = "na.or.complete") ) #si

  t1 <- mean(c(c2),  na.rm=TRUE) #sender K
  t2 <- sd(c(c2),  na.rm=TRUE) #sender K
  col1 <- abs( t1 -t2)
  t1 <-  mean(c(c1), na.rm=TRUE) #receiver K
  t2 <- sd(c(c1), na.rm=TRUE) #receiver K
  col2 <-  abs( t1 -t2)
  col3 <- mean(c(c4),  na.rm=TRUE) #sender i
  col4 <- sd(c(c4),  na.rm=TRUE) #sender i
  col5 <- mean(c(c3),  na.rm=TRUE) #receiver i
  col6 <- sd(c(c3),  na.rm=TRUE) #receiver i
  newDATA <- matrix(c(col1, col2, col3, col4, col5, col6), ncol = 6)
  out2 <- mclust::predict.MclustDA( mp, newDATA)
  best_model <- namesMod[as.numeric( out2$classification)]
  model_probabilities <- out2$z
  colnames(model_probabilities) <- namesMod

  out <- list( bestModel = best_model, modelProbs = model_probabilities,
               modClass = out_mod_class)
  class(out) <- "heuSearch"
  return( out )
}

print.heuSearch <- function(x, ...){
  cat("\n", " Selected model:",  x$bestModel, "\n")
  temp <- c( substr( x$bestModel,1,1), substr( x$bestModel,2,2))
  #first letter is the sender effect type, second letter the receiver effect type
  temp1 <- switch (temp[1],
                   N = "null",
                   C = "constant",
                   V = "variable"
  )
  temp2 <- switch (temp[2],
                   N = "null",
                   C = "constant",
                   V = "variable"
  )
  cat("\n", " Latent space model with", temp1, "sender effect and", temp2,
      "receiver effect.", "\n")
}
