#' Sequential generalized Hochberg and Hommel procedures based on group-sequential p-values
#'
#' `seqgsphh()` implements the sequential Generalized Hochberg and Hommel
#' procedures based on group-sequential p-values.
#'
#' @param pm Numeric matrix of group-sequential p-values for different
#'   hypotheses (in row) at different times (in column).
#' @param alpha Numeric scalar of the overall family-wise error rate.
#' @param epsilon Numeric scalar indicating the lower bound for `alpha`.
#' @param precision Integer scalar for precision of the values, obsolete for
#'   backward compatibility.
#' @param method "Hochberg" or "Hommel"
#' @return List with elements
#'   * `rejected`: the index set of rejected hypotheses
#'   * `decisionsm`: rejection decision for each endpoint (row) at each timepoint
#'     (column)
#'   * `cumdecisionsm`: cumulative rejection decision for each endpoint (row) at
#'     each timepoint (column)
#' @author Xiaodong Luo
#' @concept Hochberg procedure
#' @concept Hommel procedure
#' @concept group-sequential
#' @concept group-sequential p-values
#' @examples
#' pm <- matrix(rep(c(0.03, 0.04, 0.01), times = 2), ncol = 3, nrow = 2)
#' seqgsphh(pm = pm, alpha = 0.025, method = "Hochberg")
#' seqgsphh(pm = pm, alpha = 0.025, method = "Hommel")
#' @export
seqgsphh=function(pm=matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2),
                   alpha=0.025,epsilon=1.0e-10,precision=10,method='Hochberg'){
  ### Sequential Generalized Hochberg and Hommel procedures based on group-sequential p-values
  #pm: group-sequential p-values
  #alpha: FWER
  n=nrow(pm)
  s=ncol(pm)
  mseq=seq(1,n,by=1)
  decisionsm=alphaused=cumdecisionsm=matrix(0,nrow=n,ncol=s)
  
  t=1
  alphastar=rep(alpha,n)
  if (method=='Hochberg'){
    abc=hochbergd(pvalues=pm[,t],alpha=alphastar,epsilon=epsilon,precision=precision)
  }
  else if (method=='Hommel'){
    abc=hommeld(pvalues=pm[,t],alpha=alphastar,epsilon=epsilon,precision=precision)
  }
  rejected=mseq[abc$decisions]
  decisionsm[,t]=cumdecisionsm[,t]=abc$decisions
  if (s>1){
     for (t in 2:s){
       alphastar=rep(alpha,n)
       if (method=='Hochberg'){
         abc=hochbergd(pvalues=pm[,t],alpha=alphastar,epsilon=epsilon,precision=precision)
       }
       else if (method=='Hommel'){
         abc=hommeld(pvalues=pm[,t],alpha=alphastar,epsilon=epsilon,precision=precision)
       }
       rejected=union(rejected,mseq[abc$decisions])
       decisionsm[,t]=abc$decisions
       cumdecisionsm[,t]=(abc$decisions|cumdecisionsm[,t-1])
     }
  }
  if (length(rejected)>0) rejected=sort(rejected)
  list(rejected=rejected,decisionsm=decisionsm,cumdecisionsm=cumdecisionsm)
}
