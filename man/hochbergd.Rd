\name{hochbergd}
\alias{hochbergd}
\title{Hochberg procedure with different alphas for different endpoints}
\description{Hochberg procedure with different alphas for different endpoints}
\usage{
hochbergd(pvalues,alpha,epsilon=1.0e-10,precision=10)
}
\arguments{
  \item{pvalues}{p-values from different endpoints}
  \item{alpha}{same length as \code{pvalues} with (different) alphas for different endpoints}
  \item{epsilon}{lower bound for the alpha}
  \item{precision}{precision of the values}
}
\details{This procedure handles Hochberg procedure with different alphas for different endpoints.}
\value{\item{decisions}{an index of rejected hypotheses}}
\author{Xiaodong Luo}
\examples{
hochbergd(pvalues=runif(5),alpha=seq(0.01,0.025,len=5),
          epsilon=1.0e-10,precision=10)
}
\keyword{Hochberg procedure}

