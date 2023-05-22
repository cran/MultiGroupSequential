\name{spendingfun}
\alias{spendingfun}
\title{Calculate alpha spending function}
\description{This utility function calculates alpha spending. Note that the OBF and Pocock spending functions are not the originally proposed ones, they are the modified ones that closely resemble the original versions. You might still see some differences}
\usage{
spendingfun(alpha,fractions=seq(0.2,1,by=0.2),family="OBF",rho=1)
}
\arguments{
  \item{alpha}{overall alpha to be spent}
  \item{fractions}{information fractions}
  \item{family}{family of spending functions, one of "OBF", "pocock", "power"}
  \item{rho}{auxiliary parameter for OBF and power family}
}
\details{
OBF: \eqn{2\{1-\Phi(\Phi^{-1}(1-\alpha/2)/t^{\rho/2})\}};
pocock: \eqn{\alpha \log\{1+(e-1)*t\}};
power: \eqn{\alpha*t^{\rho}}
}
\value{\item{aseq}{alpha spending}
}
\author{Xiaodong Luo}
\examples{
spendingfun(alpha=0.025,fractions=seq(0.2,1,by=0.2),family="OBF",rho=1)
}
\keyword{group-sequential}
\keyword{alpha spending}

