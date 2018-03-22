\name{vickers}
\alias{vickers}
\alias{sex}
\title{Vickers-Chan 7th Graders Multidimensional Network}

\description{
The data by Vickers et al. (1981) represent three different social relations among 29 seventh grade students in a school in Victoria, Australia. A binary covariate matrix is available and it describes if two students share the same sex or not.
}

\usage{
data(vickers)
}


\details{
The data by Vickers et al. (1981) are stored in a list named \code{vickers} and represent three different social relations among 29 seventh grade students, 12 boys and 17 girls, in a school in Victoria, Australia. Each student was asked if
\itemize{
\item{\emph{Network 1}: he\ she got on with another given student;}
\item{\emph{Network 2}: he\ she considered another given student as his\ her best friend;}
\item{\emph{Network 3}: he\ she would have liked to work with another given student.}
}
The matrix \code{sex} is an edge-specific binary covariate recording the gender of pairs of students. If a pair of students is of the same sex, the corresponding entry in the matrix is 0, otherwise is 1.
}

\references{
Vickers, M. and Chan, S. (1981). Representing Classroom Social Structure. \emph{Melbourne: Victoria Institue of Secondary Education.}
}

\examples{
data(vickers)
}
