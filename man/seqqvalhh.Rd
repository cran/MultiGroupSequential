\name{seqqvalhh}
\alias{seqqvalhh}
\title{Sequential Hochberg/Hommel procedure for multiple hypotheses based on q-values}
\description{Sequential Hochberg/Hommel procedure for multiple hypotheses based on q-values}
\usage{
seqqvalhh(pm=matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2),
          alpham=matrix(rep(c(0.02,0.03,0.05),each=2),ncol=3,nrow=2),
          epsilon=1.0e-10,precision=10,method='Hochberg')          
          
}
\arguments{
  \item{pm}{a matrix of group-sequential p-values for different hypotheses (in row) at different times (in column)}
  \item{alpham}{a matrix of alpha spending corresponding to the p-values \code{pm}}
  \item{epsilon}{lower bound for the alpha}
  \item{precision}{precision of the values}
  \item{method}{"Hochberg" or "Hommel"}
}
\value{\item{rejected}{the index set of rejected hypotheses}
       \item{decisionsm}{rejection decision for each endpoint (row) at each timepoint (column)}
       \item{cumdecisionsm}{cumulative rejection decision for each endpoint (row) at each timepoint (column)}
       \item{alphaused}{alpha levels actually used for each endpoint (row) at each timepoint (column)}
       
}
\author{Xiaodong Luo}
\examples{
pm=matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2)
alpham=matrix(rep(c(0.02,0.03,0.05),each=2),ncol=3,nrow=2)
seqqvalhh(pm=pm,alpham=alpham,method="Hochberg")
seqqvalhh(pm=pm,alpham=alpham,method="Hommel")
}
\keyword{group-sequential}
\keyword{Hochberg}
\keyword{Hommel}
\keyword{q-values}


