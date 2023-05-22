\name{calgsp1}
\alias{calgsp1}
\title{Calculate group-sequential p-values for one hypothesis}
\description{This will Calculate the group-sequential p-values for one hypothesis}
\usage{
calgsp1(sx=qnorm(1-c(0.03,0.04,0.01)),scrit=qnorm(1-c(0.01,0.02,0.025)),
        salpha=c(0.01,0.02,0.025),smatrix=diag(3),sided=1)
}
\arguments{
  \item{sx}{sequence of test statistics, assumed to be multivariate normal with variance=1 and correlation matrix as \code{smatrix}}
  \item{scrit}{sequence of critical values for the statistics \code{sx}}, it should be computed beforehand
  \item{salpha}{sequence of cumulative alpha levels for the statistics \code{sx}}
  \item{smatrix}{correlation matrix of the test statistics \code{sx}}
  \item{sided}{-1: (one-sided)reject if test stat is smaller than or equal to the critical value;
                1: (one-sided)reject if test stat is greater than or equal to the critical value;
                0: (two-sided)reject if the absolute value of the test stat is greater than the critical value.}
}
\value{
  \item{pa}{group-sequential p-values}
}
\author{Xiaodong Luo}
\examples{
calgsp1(sx=qnorm(1-c(0.03,0.04,0.01)),scrit=qnorm(1-c(0.01,0.02,0.025)),
        salpha=c(0.01,0.02,0.025),smatrix=diag(3),sided=1)
}
\keyword{group-sequential p-values}



