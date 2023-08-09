##' Transform information fractions into correlation matrix
##'
##' `inftocor()` transforms information (fractions) into correlation matrix.
##'
##' @param ir Numeric vector of the sequence of information fractions. All
##'   elements should be between 0 and 1 with the last one being exactly 1.
##' @return List with an element named `cor` for the correlation matrix.
##' @author Xiaodong Luo
##' @concept correlation matrix
##' @concept information fraction
##' @examples
##' inftocor(ir = c(0.2, 0.5, 1.0))
##' @export
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

