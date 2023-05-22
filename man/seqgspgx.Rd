\name{seqgspgx}
\alias{seqgspgx}
\title{Sequential graphical procedure for multiple hypotheses based on group-sequential p-values}
\description{Sequential graphical procedure for multiple hypotheses based on group-sequential p-values}
\usage{
seqgspgx(pm=matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2),alpha=0.025,
                    W=c(0.6,0.4),G=rbind(c(0,1),c(1,0)))
}
\arguments{
  \item{pm}{a matrix of group-sequential p-values for different hypotheses (in row) at different times (in column)}
  \item{alpha}{overall familywise error rate}
  \item{W}{the weights of the graph}
  \item{G}{the transition matrix of the graph}
}
\value{\item{rejected}{the index set of rejected hypotheses}
       \item{decisionsm}{rejection decision for each endpoint (row) at each timepoint (column)}
       \item{cumdecisionsm}{cumulative rejection decision for each endpoint (row) at each timepoint (column)}
}
\author{Xiaodong Luo}
\examples{
seqgspgx(pm=matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2),alpha=0.025,
                    W=c(0.6,0.4),G=rbind(c(0,1),c(1,0)))
}
\keyword{group-sequential}

