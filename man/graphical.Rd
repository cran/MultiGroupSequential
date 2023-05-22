\name{graphical}
\alias{graphical}
\title{Graphical procedure}
\description{This function performs graphical procedure to test multiple hypotheses}
\usage{
graphical(p=c(0.01,0.04,0.03),W=c(0.5,0.25,0.25),
          G=rbind(c(0,1,0),c(0,0,1),c(1,0,0)),alpha=0.05)
}
\arguments{
  \item{p}{vector of p-values for the hypotheses}
  \item{W}{the weights of the graph}
  \item{G}{the transition matrix of the graph}
  \item{alpha}{overall type-1 error rate}
}
\value{
   \item{rej.h}{a vector indicating which hypothesis is rejected}
}
\author{Kaiyuan Hua, Xiaodong Luo}
\examples{
graphical(p=c(0.02,0.03,0.01))
}
\keyword{graphical procedure}
\keyword{multiple comparison}