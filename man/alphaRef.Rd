\name{alphaRef}
\alias{alphaRef}
\title{Compute Intercept Value in the Reference Network}

\description{
Computes the value for the intercept in the reference network.
}

\usage{
alphaRef(Y, D = 2)
}

\arguments{
\item{Y}{A three-dimensional array or list of \eqn{(n\times n)}{n x n} adjacency matrices composing the multidimensional network. A list will be converted to an array. If an array, the dimension of \code{Y} must be \code{(n,n,K)}, where \code{n} is the number of nodes and \code{K} the number of networks. Structural missing values (\code{NA}) are allowed, but no imputation is performed; see details.}

\item{D}{The dimension of the latent space, with \code{D > 0}. The default value is \code{D = 2}.}

}

\value{
The function returns the intercept value in the reference network.
}

\details{
The function computes the value for the intercept in the reference network (first network) of the multiplex. It is calculated taking into account the average effect of the latent space on edge probabilities (approximated with the constant \eqn{2}{2}) and the observed mean probability of an edge in the first network (\eqn{Y^{(1)}}{Y_1}):
\deqn{p^{(1)}=\frac{\sum_{i=1}^n \sum_{j =1}^n y_{ij}^{(1)}}{n (n-1)},}{p_1 = sum( Y_1) \ ( n * ( n - 1 ) ),}
with \eqn{n}{n} the number of nodes in the network.
Then, the reference intercept is computed as:
\deqn{\log \Bigl( \frac{p^{(1)}}{1 -p^{(1)}} \Bigr) + 2.}{ log ( p_1 \ ( 1 - p_1) ) + 2.}
}


\references{
D'Angelo, S. and Murphy, T. B. and Alfò, M. (2018). Latent space modeling of multidimensional networks with application to the exchange of votes in the Eurovision Song Contest. \href{https://arxiv.org/abs/1803.07166}{arXiv}.
}

\examples{
data(vickers)
alphaR <- alphaRef(vickers, D = 2)  # model without node-specific effect
}

\seealso{
\code{\link{multiNet}}
}

