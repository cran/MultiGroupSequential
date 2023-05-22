\name{seqgsphh}
\alias{seqgsphh}
\title{Sequential Hochberg/Hommel procedure for multiple hypotheses based on group-sequential p-values}
\description{Sequential Hochberg/Hommel procedure for multiple hypotheses based on group-sequential p-values}
\usage{
seqgsphh(pm=matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2),
                   alpha=0.025,epsilon=1.0e-10,precision=10,method='Hochberg')

}
\arguments{
  \item{pm}{a matrix of group-sequential p-values for different hypotheses (in row) at different times (in column)}
  \item{alpha}{overall familywise error rate}
  \item{epsilon}{lower bound for the alpha}
  \item{precision}{precision of the values}
  \item{method}{"Hochberg" or "Hommel"}
}
\value{\item{rejected}{the index set of rejected hypotheses}
       \item{decisionsm}{rejection decision for each endpoint (row) at each timepoint (column)}
       \item{cumdecisionsm}{cumulative rejection decision for each endpoint (row) at each timepoint (column)}
}
\author{Xiaodong Luo}
\examples{
pm=matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2)
seqgsphh(pm=pm,alpha=0.025,method="Hochberg")
seqgsphh(pm=pm,alpha=0.025,method="Hommel")
}
\keyword{group-sequential}
\keyword{Hochberg}
\keyword{Hommel}
\keyword{group-sequential p-values}

