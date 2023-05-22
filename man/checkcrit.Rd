\name{checkcrit}
\alias{checkcrit}
\title{Check critical values}
\description{This utility function checks the critical values}
\usage{
checkcrit(scrit=qnorm(c(0.01,0.02,0.025)),salpha=c(0.01,0.02,0.025),
          smatrix=diag(3),sided=1)
}
\arguments{
  \item{scrit}{sequence of critical values}
  \item{salpha}{sequence of cumulative alpha levels}
  \item{smatrix}{general correlation matrix}
  \item{sided}{-1: (one-sided)reject if test stat is smaller than or equal to the critical value;
                    1: (one-sided)reject if test stat is greater than or equal to the critical value;
                    0: (two-sided)reject if the absolute value of the test stat is greater than the critical value}
}
\value{
    \item{crit.value}{critical values}
    \item{salpha}{sequence of cumulative alpha levels same as input}
}
\author{Xiaodong Luo}
\examples{
checkcrit(scrit=qnorm(c(0.01,0.02,0.025)),salpha=c(0.01,0.02,0.025),
          smatrix=diag(3),sided=1)
}
\keyword{group-sequential}
\keyword{critical values}
\keyword{efficacy boundary}

