\name{multiNet}
\alias{multiNet}
\alias{print.multiNet}
\title{Latent Space Models for Multivariate Networks}

\description{
Implements latent space models for multivariate networks (multiplex) via MCMC algorithm.
}

\usage{
multiNet(Y, niter = 1000, D = 2,
         muA = 0, tauA = NULL, nuA = 3,
         muB = 0, tauB = NULL, nuB = 3,
         muL = 0, tauL = NULL, nuL = 3,
         alphaRef = 0.1,
         covariates = NULL,
         DIC = FALSE,
         burnIn = round(niter*0.3),
         trace = TRUE,
         allChains = FALSE,
         refSpace = NULL)
}


\arguments{
\item{Y}{A three-dimensional array or list of \eqn{(n\times n)}{n x n} adjacency matrices composing the multidimensional network. A list will be converted to an array. If an array, the dimension of \code{Y} must be \code{(n,n,K)}, where \code{n} is the number of nodes and \code{K} the number of networks. Missing values (\code{NA}) are allowed; see details.}

\item{niter}{The number of MCMC iterations. The default value is \code{niter = 1000}.}

\item{D}{The dimension of the latent space, with \code{D > 0}. The default value is \code{D = 2}.}

\item{muA, muB, muL}{Mean hyperparameters, see details.}

\item{tauA, tauB, tauL}{Mean hyperparameters, see details.}

\item{nuA, nuB, nuL}{Variance hyperparameters, see details.}

\item{alphaRef}{The value for the intercept in the first network (the reference network). This value can be specified by the user on the basis of prior knowledge or computed using the function \code{\link{alphaRef}}, see details.}

\item{covariates}{An array or a list with edge-covariates matrices. A list is automatically converted to an array. Covariates can be either continuous or discrete and must be constant throughout the views of the multiplex. The dimension of \code{covariates} is \code{(n,n,L)}, where \code{n} is the number of nodes and \code{L} the number of covariates, that is, the number of covariates matrices. Missing values (\code{NA}) are not allowed.}

\item{DIC}{A logical value indicating wether the DIC (Deviance Information Criterion) should be computed. The default is \code{DIC = FALSE}.}

\item{burnIn}{A numerical value, the number of iterations of the chain to be discarded when computing the posterior estimates. The default value is \code{burnIn = round(niter*0.3)}.}

\item{trace}{A logical value indicating if a progress bar should be printed.}

\item{allChains}{A logical value indicating if the full parameter chains should also be returned in output. The default value is \code{allChains = FALSE}.}

\item{refSpace}{Optional. A matrix containing a set of reference values for the latent coordinates of the nodes. Its dimension must be \code{(n, D)}, where \code{n} is the number of nodes and \code{D} the number of dimensions of the latent space. The coordinates stored in the matrix \code{refSpace} are compared with the estimated ones at each iteration via Procrustes correlation. High values of the correlation index indicate that the estimated coordinates are a translation and/or a rotation of the coordinates in \code{refSpace}.}
}


\value{
An object of class \code{'multiNet'} containing the following components:
\item{n}{The number of nodes in the multidimensional network.}
\item{K}{The number of networks in the multidimensional network.}
\item{D}{The number of dimensions of the estimated latent space.}
\item{parameters}{A list with the following components:
\itemize{
\item{\code{alpha}}{ is a list with two components: the means of the posterior distributions and the standard deviations of the posterior distributions for the intercept parameters; }
\item{\code{beta}}{ is a list with two components: the means of the posterior distributions and the standard deviations of the posterior distributions for the latent space coefficient parameters;}
\item{\code{lambda}}{ is a list with two components: the means of the posterior distributions and the standard deviations of the posterior distributions for the covariate coefficient parameters.}
}
}
\item{latPos}{A list with posterior estimates of means and standard deviations of the latent coordinates.}
\item{accRates}{A list with the following components:
\itemize{
\item{\code{alpha}}{ is a vector with the acceptance rates for the intercept parameters;}
\item{\code{beta}}{ is a vector with the acceptance rates for the latent space coefficient parameters;}
\item{\code{lambda}}{ is a vector with the acceptance rates for the covariate coefficient parameters;}
\item{\code{latPos}}{ is a vector with the acceptance rates for the latent coordinates of the nodes.}
}
}
\item{DIC}{The Deviance Information Criterion of the estimated model. Computed only if \code{DIC = TRUE} in input.}
\item{allChains}{If \code{allChains = TRUE}, a list with the following components is returned:
\itemize{
\item{\code{parameters}}{ is a list with the estimated posterior distributions of the model parameters: \eqn{\alpha}{\alpha}, \eqn{\beta}{\beta} and \eqn{\lambda}{\lambda};}
\item{\code{latPos}}{ is an array with the posterior distributions of the latent coordinates of each node;}
\item{\code{priorParameters}}{ is a list with the estimated posterior distributions of the parameters of the prior distributions of \eqn{\alpha}{\alpha}, \eqn{\beta}{\beta} and \eqn{\lambda}{\lambda}. }
}
}
\item{corrRefSpace}{A numerical vector containing the values of the Procrustes correlation between the reference space and the estimated one, computed at each mcmc iteration. Only outputed when \code{refSpace} is given, otherwise \code{NULL}.}
\item{info}{A list with some information on the estimated model:
\itemize{
\item{\code{call}}{ contains the function call;}
\item{\code{niter}}{ is the number of MCMC iterations;}
\item{\code{burnIn}}{ is the number of initial iterations to discarded when computing the estimates;}
\item{\code{covariates}}{ is the covariates array, if present;}
\item{\code{L}}{ is the number of covariates.}
}
}
}

\details{
The function estimates a latent space model for multidimensional networks (multiplex) via MCMC. The model assumes that the probability of observing an arc between any two nodes is inversely related to their distance in a low-dimensional latent space. Hence, nodes close in the latent space have a higher probability of being connected across the views of the multiplex than nodes far apart.

The probability of an edge beteween nodes \eqn{i}{i} and \eqn{j}{j} in the \eqn{k^{th}}{k_{th}} network is defined as:

\deqn{
 P \Bigl( y_{ij}^{(k)} = 1 | \Omega^{(k)} , d_{ij}, \lambda \Bigr)=
\frac{ C_{ij}^{(k)} }{1 +  C_{ij}^{(k)} }.}{
P ( y_{ijk} = 1 | \Omega_k , d_{ij} , \lambda ) =  C_{ijk} \ ( 1 + C_{ijk} ).}

with \eqn{C_{ij}^{(k)} = \exp \{\alpha^{(k)}-\beta^{(k)} d_{ij} -\lambda x_{ij} \} }{C_{ijk} = exp( \alpha_k - \beta_k * d_{ij} - \lambda * x_{ij} )}. The arguments of \eqn{C_{ij}^{(k)}}{C_{ijk}} are:
\itemize{
\item{}{The squared Euclidean distance between nodes \eqn{i}{i} and \eqn{j}{j} in the latent space, \eqn{d_{ij}}{d_{ij}.}}
\item{}{A coefficient \eqn{\lambda}{\lambda} to scale the edge-specific covariate \eqn{x_{ij}}{x_{ij}}. If more than one covariate is introduced in the model, their sum is considered, with each covariate being rescaled by a specific coefficient \eqn{\lambda_l}{\lambda_l}. Edge-specific covariates are assumed to be inversely related to edge probabilities, hence \eqn{\lambda \geq 0}{ \lambda => 0 }.}
\item{}{A vector of network-specific parameters, \eqn{\Omega^{(k)} = (\alpha^{(k)},\beta^{(k)}) }{\Omega_k = ( \alpha_k, \beta_k )}. These parameters are:
\itemize{
\item{}{A rescaling coefficient \eqn{\beta^{(k)} }{\beta_k}, which weights the importance of the latent space in the \eqn{k^{th}}{k_{th}} network, with \eqn{\beta^{(k)} \geq 0 }{\beta_k => 0}. In the first network (that is the reference network), the coefficient is fixed to \eqn{\beta^{(1)} = 1}{\beta_1 = 1} for identifiability reasons.  }
\item{}{An intercept parameter \eqn{\alpha^{(k)} }{\alpha_k}, which corresponds to the largest edge probability allowed in the \eqn{k^{th}}{k_{th}} network. Indeed, when  \eqn{\beta^{(k)} = 0 }{\beta_k = 0} and when no covariate is included, the probability of having a link between a couple of nodes is that of the random graph:
\deqn{
 P \Bigl( y_{ij}^{(k)} = 1 |  \alpha^{(k)} \Bigr)=
\frac{ \exp \{ \alpha^{(k)}\} }{1 +  \exp \{\alpha^{(k)}\} }.}{
P ( y_{ijk} = 1 | \alpha_k ) = exp( \alpha_k ) \ ( 1 + exp( \alpha_k ) ).}
}
The intercepts have a lower bound corresponding to \eqn{\log \Bigl( \frac{\log (n)}{ n - \log(n)} \Bigr) }{ log ( log( n ) \ ( n - log( n )  ) )}. For identifiability reasons, the intercept of the first network needs to be fixed. Its value can be either specified by the user on the basis of prior knowledge or computed with the function \code{\link{alphaRef}}.
}
}
}
Inference on the model parameters is carried out via a MCMC algorithm. A hierarchical framework is adopted for estimation, where the parameters of the distributions of \eqn{\alpha}{\alpha}, \eqn{\beta}{\beta} and \eqn{\lambda}{\lambda} are considered nuisance parameters and assumed to follow hyper-prior distributions. The parameters of these hyperpriors need to be fixed and are the following:
\itemize{
\item{}{\code{tauA, tauB} and \code{tauL} are the scale factors for the variances of the hyperprior distributions for the mean parameters of \eqn{\alpha^{(k)}, \beta^{(k)}}{\alpha_k , \beta_k} and \eqn{\lambda_l}{\lambda_l}. If not specified by the user, \code{tauA} and \code{tauB} are computed as \eqn{(K-1)\ K }{( K - 1 ) \ K }, if \eqn{K > 1}{K > 1}, otherwise they are set to \eqn{0.5}{0.5}. Parameter \code{tauL} is calculated as \eqn{(L-1)\ K }{( L - 1 ) \ K }, if \eqn{L > 1}{L > 1}, otherwise it is set to \eqn{0.5}{0.5}.    }
\item{}{\code{muA, muB} and \code{muL} are the means of the hyperprior distributions for the mean parameters of \eqn{\alpha^{(k)}, \beta^{(k)}}{\alpha_k , \beta_k} and \eqn{\lambda_l}{\lambda_l}. If not specified by the user, they are all set to \eqn{0}{0}. }
\item{}{\code{nuA, nuB} and \code{nuL} are the degrees of freedom of the hyperprior distributions for the variance parameters of \eqn{\alpha^{(k)}, \beta^{(k)}}{\alpha_k , \beta_k} and \eqn{\lambda_l}{\lambda_l}. If not specified by the user, they are all set to \eqn{3}{3}. }
}

Missing data are considered structural and correspond to edges missing because one or more nodes are not observable in some of the networks of the multiplex. No imputation is performed, instead, the term corresponding to the missing edge is discarded in the computation of the likelihood function. For example, if either node \eqn{(i)}{(i)} or \eqn{(j)}{(j)} is not observable in network \eqn{(k)}{(k)}, the edge \eqn{(i,j)}{(i,j)} is missing and the likelihood function for network \eqn{(k)}{(k)} is calculated discarding the corresponding \eqn{(i,j)}{(i,j)} term. Notice that the model assumes a single common generative latent space for the whole multidimensional network. Thus, discarding the \eqn{(i,j)}{(i,j)} term in the \eqn{k^{th}}{k_{th}} network does not prevent from recovering the coordinates of nodes \eqn{i}{i} and \eqn{j}{j} in the latent space.
}

%\author{Silvia D'Angelo, Michael Fop}


\examples{
data(vickers)
alphaR <- alphaRef(vickers, D = 2)  # compute alpha reference value

it <- 10     # small number of iterations just for example

# 2-dimensional latent space model, no covariates
mod <- multiNet(vickers, niter = it, D = 2,
                alphaRef = alphaR)

# 2-dimensional latent space model, sex as covariate
mod <- multiNet(vickers, niter = it, D = 2,
                alphaRef = alphaR,
                covariates = sex)
}

\references{
D'Angelo, S. and Murphy, T. B. and Alfò, M. (2018). Latent space modeling of multidimensional networks with application to the exchange of votes in the Eurovision
Song Contest. \href{https://arxiv.org/abs/1803.07166}{arXiv}.

}

\seealso{
\code{\link{alphaRef}}
}
