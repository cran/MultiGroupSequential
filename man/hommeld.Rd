\name{hommeld}
\alias{hommeld}
\title{Hommel procedure with different alphas for different endpoints}
\description{Hommel procedure with different alphas for different endpoints}
\usage{
hommeld(pvalues,alpha,epsilon=1.0e-10,precision=10)
}
\arguments{
  \item{pvalues}{p-values from different endpoints}
  \item{alpha}{same length as \code{pvalues} with (different) alphas for different endpoints}
  \item{epsilon}{lower bound for the alpha}
  \item{precision}{precision of the values, obsolete for backward compatibility}
}
\details{It turns out [hommel](https://cran.r-project.org/web/packages/hommel/) can handle Hommel procedure with different alpha's for different endpoints, the function 'hommeld' is just a wrapper function. }
\value{\item{decisions}{an index of rejected hypotheses}}
\author{Xiaodong Luo}
\examples{
hommeld(pvalues=runif(5),alpha=seq(0.01,0.025,len=5),
        epsilon=1.0e-10,precision=10)
}
\keyword{Hommel procedure}

