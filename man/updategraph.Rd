\name{updategraph}
\alias{updategraph}
\title{Update graph}
\description{This function updates the graph when only a subset of original hypotheses is concerned}
\usage{
updategraph(S1=c(2,3),W0=c(0.5,0.5,0,0),
            G0=rbind(c(0,0,1,0),c(0,0,0,1),c(0,1,0,0),c(1,0,0,0)),
            S0=seq(1,length(W0),by=1))
}
\arguments{
  \item{S1}{the subset of hypotheses, S1 must be a non-empty subset of \code{S0} and must be sorted increasingly}
  \item{W0}{the weights of the graph}
  \item{G0}{the transition matrix of the graph}
  \item{S0}{the set of hypotheses from 1 to n}
}
\value{
   \item{S1}{same as the input \code{S1}}
   \item{W1}{the weights of the updated graph}
   \item{G1}{the transition matrix of the updated graph}
}
\author{Xiaodong Luo}
\examples{
## We can use the above function to produce a closed testing tree
# A function to create power set
powerset <- function(x) {
  sets <- lapply(1:(length(x)), function(i) combn(x, i, simplify = FALSE))
  unlist(sets, recursive = FALSE)
}

n=3 #number of hypotheses
pn=2^n-1
pset=powerset(seq(1,n,by=1)) #create the power set
df <- data.frame(matrix(ncol =1+n, nrow = 0)) #create the dataset
colnames(df) <- c("Test",paste0("H",seq(1,n,by=1),sep=""))

W0=c(1/3,1/3,1/3) #the weights of the graph 
m=rbind(H1=c(0,1/2,1/2),
        H2=c(1/2,0,1/2),
        H3=c(1/2,1/2,0))
G0=matrix(m,nrow=3,ncol=3) #the transition matrix of the graph

for (j in 1:pn){
    abc=updategraph(S1=pset[[j]],W0=W0,G0=G0)
    temp=rep("-",n)
    temp[pset[[j]]]=abc$W1
    temp=c(paste(pset[[j]],collapse=""),temp)
    df[j,]=temp
}
df #the dataframe lists the closed testing tree
}
\keyword{graphical procedure}
\keyword{multiple comparison}
\keyword{closed testing}


