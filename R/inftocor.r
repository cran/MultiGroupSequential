#Transform information (fractions) into correlation matrix
inftocor=function(ir=c(0.2,0.5,1)){
  # ir: sequence of information or information fractions
  nir=length(ir)
  dcorr=diag(nir)
  
  for (u in 1:nir){
    for (v in u:nir){
      dcorr[u,v]=dcorr[v,u]=sqrt(ir[u]/ir[v])
    }
  }
  list(cor=dcorr)
}

