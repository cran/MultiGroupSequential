\name{seqmbgx}
\alias{seqmbgx}
\title{Maurer-Bretz sequential graphical approach}
\description{This will conduct group-sequential testing for multiple hypotheses based on Maurer-Bretz approach}
\usage{
seqmbgx(xm=qnorm(matrix(rep(c(0.03,0.04,0.01),times=4),ncol=3,nrow=4)),
        informationm=matrix(rep(c(0.4,0.8,1),each=4),ncol=3,nrow=4),
        spending=rep("OBF",4),param.spending=rep(1,4),
        alpha=0.025,sided=-1,
        W=c(0.5,0.5,0,0),G=rbind(c(0,0,1,0),c(0,0,0,1),c(0,1,0,0),c(1,0,0,0)),
        tol=1e-10,retrospective=0)
}
\arguments{
  \item{xm}{a matrix of test statistics for each endpoint (in row) and each interim (in column)}
  \item{informationm}{a matrix of information fractions for the statistics \code{xm}}
  \item{spending}{spending functions for each endpoint}
  \item{param.spending}{parameters in the spending functions}
  \item{alpha}{overall familywise error rate}
  \item{sided}{    -1: (one-sided)reject if test stat is smaller than or equal to the critical value;
                    1: (one-sided)reject if test stat is greater than or equal to the critical value;
                    0: (two-sided)reject if the absolute value of the test stat is greater than the critical value}
  \item{W}{the weights of the graph}
  \item{G}{the transition matrix of the graph}
  \item{tol}{tolerance level for computing the critical values}
  \item{retrospective}{retrospective: 0 (default) only compares the current test statistic with the updated critical value, 1 compares all the test statistics up to the current one with the updated critical values. Even though retrospective looking at the values is statistically valid in terms of control of the type-1 error rate, no retrospective looking at the past comparisons avoids the dilemma of retrospectively increasing the alpha level for the un-rejected hypothesis in the past}
}
\value{\item{Hrej}{rejected hypotheses}
       \item{rejected}{the index set of rejected hypotheses}
       \item{decisionsm}{rejection decision for each endpoint (row) at each timepoint (column)}
       \item{cumdecisionsm}{cumulative rejection decision for each endpoint (row) at each timepoint (column)}
}
\author{Xiaodong Luo}
\examples{
xm=qnorm(matrix(rep(c(0.03,0.04,0.01),times=4),ncol=3,nrow=4))
im=matrix(rep(c(0.4,0.8,1),each=4),ncol=3,nrow=4)
spending=rep("OBF",4)
param.spending=rep(1,4)
W=c(0.5,0.5,0,0)
G=rbind(c(0,0,1,0),c(0,0,0,1),c(0,1,0,0),c(1,0,0,0))
seqmbgx(xm=xm,informationm=im,spending=spending,param.spending=param.spending,
          alpha=0.025,W=W,G=G,retrospective=0)
}
\keyword{graphical procedure}
\keyword{group-sequential}
\keyword{Maurer-Bretz}

