\name{spaceNet}
\alias{spaceNet}
\alias{spaceNet-package}
\title{Latent Space Models for Multivariate Networks}

\description{
A package for latent space models for binary multivariate networks (multiplex).
The model assumes that the nodes in the multiplex lie in a low-dimensional latent space. The probability of two nodes being connected is inversely related to their distance in this latent space: nodes close in the space are more likely to be linked, while nodes that are far apart are less likely to be connected. The model allows the inclusion of node-specific sender and receiver effects and edge-specific covariates. Inference is carried out via a MCMC algorithm.}


\details{ %%%to do
The main function is \code{\link{multiNet}}, which estimates the latent space model via MCMC algorithm. Data can be inputed either as a list or an array. Also, edge-specific covariates in the form of a list or an array can be included in the model.
}

\author{Silvia D'Angelo and Michael Fop.

Mantainer: Silvia D'Angelo \email{silvia.dangelo@ucd.ie}

}


\references{
D'Angelo, S. and Murphy, T. B. and Alfò, M. (2018). Latent space modeling of multidimensional networks with application to the exchange of votes in the Eurovision
Song Contest. \href{https://arxiv.org/abs/1803.07166}{arXiv}.


D'Angelo, S. and Alfò, M. and Murphy, T. B. (2018). Node-specific effects in latent space modelling of multidimensional networks. \href{https://arxiv.org/abs/1807.03874}{arXiv}.
}
