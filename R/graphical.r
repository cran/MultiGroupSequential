### Graphical procedure
graphical=function(p=c(0.01,0.04,0.03),W=c(0.5,0.25,0.25),G=rbind(c(0,1,0),c(0,0,1),c(1,0,0)),alpha=0.05){
  #p vector of p-values
  #W weights of length(p) of the graph
  #G transition matrix of the graph
  #alpha overall type-1 error rate
  np=length(p)
  I=seq(1,np,by=1)
  have.rej <- TRUE
  ps=p
  while (have.rej){
    rej.judge <- (ps <= W * alpha)
    have.rej <- sum(rej.judge) > 0
    if (have.rej){
      rej.pos <- which(rej.judge == TRUE)
      j <- rej.pos[1]  # pick the first rejected one
      W <- W + W[j] * G[j,]  # update W
      G.new <- G
      for (l in 1:nrow(G.new)){  # update G
        for (k in 1:ncol(G.new)){
          G.new[l,k] <- (G[l,k] + G[l,j]*G[j,k])/(1 - G[l,j]*G[j,l])*(1 - (l == k))
        }
      }
      # update I <- I/{j}
      I <- I[-j]
      W <- W[-j]
      if (length(I) > 0){
        ps <- ps[-j]
        G <- as.matrix(G.new[-j,-j])
      }
      else {have.rej <- FALSE}
    }
  }
  rej.h <- rep(1,length(p))
  rej.h[I] <- 0
  # return a vector of length of hypotheses number
  # 1: rejected, 0: not rejected
  list(rej.h=rej.h)
}