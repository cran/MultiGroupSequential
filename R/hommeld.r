
hommeld=function(pvalues,alpha,epsilon=1.0e-10,precision=10){
  ## Hommel procedure that can handle different alphas for different endpoints
  alpha1=pmax(alpha,epsilon)
  qvalues=pvalues/alpha1
  #qvalues=round(qvalues,digits=precision)
  hom2 <- hommel::hommel(qvalues, simes = TRUE)
  decisions=(hommel::p.adjust(hom2) <= 1)
  list(decisions=decisions)
}