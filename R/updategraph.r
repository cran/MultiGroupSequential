### Update graph: this function is to produce the closed testing tree
updategraph=function(S1=c(2,3),W0=c(0.5,0.5,0,0),G0=rbind(c(0,0,1,0),c(0,0,0,1),c(0,1,0,0),c(1,0,0,0)),S0=seq(1,length(W0),by=1)){
  #S1: new set of hypotheses #S1 must be a non-empty subset of S0, S1 must be sorted increasingly
  #W0: initial weights
  #G0: initial matrix of length(W0)*length(W0)
  #S0: initial set of hypotheses from 1 to n
  SS=setdiff(S0,S1) 
  nss=length(SS)
  WT=W0;GT=G0
  if (nss>0){
    for (j in 1:nss){
      sj=SS[j]+1-j # this ensures we select the correct index after the graph is updated 
      WT=WT + WT[sj] * GT[sj,]  # update W 
      GT.new=GT
      for (l in 1:nrow(GT.new)){  # update G
        for (k in 1:ncol(GT.new)){
          GT.new[l,k] <- (GT[l,k] + GT[l,sj]*GT[sj,k])/(1 - GT[l,sj]*GT[sj,l])*(1 - (l == k))
        }
      }
      WT <- WT[-sj]
      GT <- as.matrix(GT.new[-sj,-sj])
    }
  }
  list(S1=S1,W1=WT,G1=GT)
}

