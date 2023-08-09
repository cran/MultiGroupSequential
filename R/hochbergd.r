##' Hochberg procedure
##'
##' `hochbergd()` computes the Hochberg procedure with different alphas
##' for different endpoints.
##'
## TODO Add more details about the method
##' @param pvalues Numeric vector of p-values from different endpoints.
##' @param alpha Numeric vector of alpha values for the different endpoints.
##'   Vector must be same length as `pvalues`.
##' @param epsilon Numeric scalar indicating the lower bound for alpha.
## TODO Change `precision` to `digits`? Or even remove this argument altogether?
##' @param precision Integer scalar of the desired number of digits to be used.
## TODO Only `decisions` element in the returned list is documented.
## What about `sqvalue` and `sq`?
## TODO Explain the returned `decision` vector further
## TODO Add info about `sqvalue` and `sq` elements in the returned list
##' @return List with element named `decisions` containing an index of rejected
##'   hypotheses.
##' @author Xiaodong Luo
##' @concept Hochberg procedure
##' @examples
##' hochbergd(
##'   pvalues = runif(5),
##'   alpha = seq(0.01, 0.025, len = 5),
##'   epsilon = 1.0e-10,
##'   precision = 10
##' )
##' @export
hochbergd=function(pvalues,alpha,epsilon=1.0e-10,precision=10){
  # Hochberg procedure that can handle different alpha's for different endpoints
  m=length(pvalues);mseq=seq(1,m,by=1)
  alpha1=pmax(alpha,epsilon)
  qvalues=pvalues/alpha1
  qvalues=round(qvalues,digits=precision)
  sq=sort(qvalues,decreasing=TRUE)
  seqc=1/(mseq)
  seqc=round(seqc,digits=precision)
  ax=(sq<=seqc)
  if (sum(ax)==0){decisions=rep(FALSE,m);sqvalue=0}
  else {istar=min(mseq[ax]);sqvalue=sq[istar];decisions=qvalues<=sqvalue}
  #list(decisions=decisions,sqvalue=sqvalue,sq=sq,seqc=seqc)
  list(decisions=decisions,sqvalue=sqvalue,sq=sq)
}
