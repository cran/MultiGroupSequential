\name{calgspn}
\alias{calgspn}
\title{Calculate group-sequential p-values for multiple hypotheses}
\description{This will Calculate the group-sequential p-values for multiple hypotheses}
\usage{
calgspn(xm=qnorm(matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2)),
        alpham=matrix(rep(c(0.02,0.03,0.05),each=2),ncol=3,nrow=2),
        critm=matrix(rep(qnorm(c(0.02,0.03,0.05)),each=2),ncol=3,nrow=2),
        matrix.list=list(diag(3),diag(3)),sided=rep(-1,2))
}
\arguments{
  \item{xm}{a matrix of test statistics for hypothesis (in row) and each interim (in column)}
  \item{alpham}{a matrix of cumulative alpha levels for the statistics \code{xm}}
  \item{critm}{a matrix of critical values for the statistics \code{xm}}, it should be computed beforehand
  \item{matrix.list}{list of correlation matrix corresponding to each hypothesis}
  \item{sided}{A vector indicating the direction of each test such that    
                   -1: (one-sided)reject if test stat is smaller than or equal to the critical value;
                    1: (one-sided)reject if test stat is greater than or equal to the critical value;
                    0: (two-sided)reject if the absolute value of the test stat is greater than the critical value.}
}
\value{
  \item{pm}{group-sequential p-values}
}
\author{Xiaodong Luo}
\examples{
calgspn(xm=qnorm(matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2)),
        alpham=matrix(rep(c(0.02,0.03,0.05),each=2),ncol=3,nrow=2),
        critm=matrix(rep(qnorm(c(0.02,0.03,0.05)),each=2),ncol=3,nrow=2),
        matrix.list=list(diag(3),diag(3)),sided=rep(-1,2))
}
\keyword{group-sequential p-values}



