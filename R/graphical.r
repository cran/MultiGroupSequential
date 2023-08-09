##' Graphical procedure
##'
##' `graphical()` performs graphical procedure to test multiple hypotheses
##'
## TODO Add more details
##' @param p Numeric vector of p-values for the hypotheses.
##' @param W Numeric vector of weigths of the graph. Must have the same length
##'   as `p`.
##' @param G Matrix of the transition matrix of the graph.
##' @param alpha Numeric scalar with the overall type-1 error rate.
## TODO Return the vector itself instead of nested within a list
##' @return A list with a single element containing a vector indicating whether
##'   hypotheses are rejected (`1`) or not (`0`).
##' @author Kaiyuan Hua, Xiaodong Luo
##' @concept graphical procedure
##' @concept multiple comparison
##' @examples
##' graphical(p = c(0.02, 0.03, 0.01))
##' @export
## TODO Change to a more informative name
graphical=function (p=c(0.01,0.04,0.03),W=c(0.5,0.25,0.25),G=rbind(c(0,1,0),c(0,0,1),c(1,0,0)),alpha = 0.05){
  np = length(p)
  I = seq(1, np, by = 1)
  have.rej <- TRUE
  ps = p
  WT=W
  GT=G
  while (have.rej) {
    rej.judge <- (ps <= WT * alpha)
    have.rej <- sum(rej.judge) > 0
    if (have.rej) {
      rej.pos <- which(rej.judge == TRUE)
      j <- rej.pos[1] # pick the first rejected one
      WT <- WT + WT[j] * GT[j, ] # update W
      GT.new <- GT
      for (l in 1:nrow(GT.new)) { # update G
        temp=(GT[l,j] * GT[j,l])
        for (k in 1:ncol(GT.new)) {
          if (l==k | temp>=1){GT.new[l,k]=0}
          else {GT.new[l, k]=(GT[l, k] + GT[l, j]*GT[j, k])/(1 - temp)}
        }
      }
      # update I <- I/{j}
      I <- I[-j]
      WT <- WT[-j]
      if (length(I) > 0) {
        ps <- ps[-j]
        GT <- as.matrix(GT.new[-j, -j])
      }
      else {
        have.rej <- FALSE
      }
    }
  }
  rej.h <- rep(1, length(p))
  rej.h[I] <- 0
  # return a vector of length of hypotheses number
  # 1: rejected, 0: not rejected
  list(rej.h = rej.h)
}



