\name{crosslist}
\alias{crosslist}
\title{Cross list all the scenarios}
\description{Cross list all the scenarios}
\usage{
crosslist(b=list(a1=c(2,3),a2=c(2,4),a3=c(0,1)))
}
\arguments{
  \item{b}{A list of lists}
}
\value{\item{df}{dataframe consisting of all the scenarios row by row}}
\author{Xiaodong Luo}
\examples{
crosslist(b=list(pbo.hazard=c(0.1,0.15,0.3), 
                 hazard.ratio=c(0.6,0.7,0.8), 
                 censoring.rate=c(0.05,0.06,0.07)
                 )
         )
}
\keyword{cross listing}

