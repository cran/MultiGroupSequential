\name{findcrit}
\alias{findcrit}
\title{Calculate critical values}
\description{This utility function calculates the critical values}
\usage{
findcrit(salpha=c(0.01,0.02,0.025),smatrix=diag(3),
         sided=1,tol=1e-10,alpha.tol=1e-11)
}
\arguments{
  \item{salpha}{sequence of cumulative alpha levels}
  \item{smatrix}{general correlation matrix}
  \item{sided}{-1: (one-sided)reject if test stat is smaller than or equal to the critical value;
                    1: (one-sided)reject if test stat is greater than or equal to the critical value;
                    0: (two-sided)reject if the absolute value of the test stat is greater than the critical value}
  \item{tol}{tolerance level for computing the critical values}
  \item{alpha.tol}{if the alpha increment is less than this, then the critical value is set to a large number determined by alpha.tol}
}
\value{\item{crit.value}{critical values}
}
\author{Xiaodong Luo}
\examples{
findcrit(salpha=c(0.01,0.02,0.025),smatrix=diag(3),
         sided=1,tol=1e-10,alpha.tol=1e-11)
}
\keyword{group-sequential}
\keyword{critical values}
\keyword{efficacy boundary}

