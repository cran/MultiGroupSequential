#' Sequential generalized Hochberg and Hommel procedures based on q-values
#'
#' @param pm Matrix of group-sequential p-values for different hypotheses (in row)
#'   at different times (in column).
#' @param alpham Matrix of alpha spending corresponding to the p-values \code{pm}.
#'   For each row, alpha levels must be non-decreasing.
#' @param epsilon Numeric scalar indicating the lower bound for alpha.
#' @param precision Integer scalar for precision of the values, obsolete for
#'   backward compatibility.
#' @param method Character scalar "Hochberg" or "Hommel".
#' @return List with elements
#'   - `rejected`: the index set of rejected hypotheses
#'   - `decisionsm`: rejection decision for each endpoint (row) at each
#'     timepoint (column)
#'   - `cumdecisionsm`: cumulative rejection decision for each endpoint (row) at
#'     each timepoint (column);
#'   - `alphaused`: alpha levels actually used for each endpoint (row) at each
#'     timepoint (column).
#' @author Xiaodong Luo
#' @concept Hochberg procedure
#' @concept Hommel procedure
#' @concept group-sequential
#' @concept q-values
#' @examples
#' pm <- matrix(rep(c(0.03, 0.04, 0.01), times = 2), ncol = 3, nrow = 2)
#' alpham <- matrix(rep(c(0.02, 0.03, 0.05), each = 2), ncol = 3, nrow = 2)
#' seqqvalhh(pm = pm, alpham = alpham, method = "Hochberg")
#' seqqvalhh(pm = pm, alpham = alpham, method = "Hommel")
#' @export
seqqvalhh=function(pm=matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2),
                 alpham=matrix(rep(c(0.02,0.03,0.05),each=2),ncol=3,nrow=2),
                 epsilon=1.0e-10,precision=10,method='Hochberg'){
  #alpham: matrix of alpha-levels, 
  n=nrow(pm)
  s=ncol(pm)
  mseq=seq(1,n,by=1)
  decisionsm=alphaused=cumdecisionsm=matrix(0,nrow=n,ncol=s)
  t=1
  alphastar=rep(1,n)
  alphaused[,t]=alphastar=alpham[,t]
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
       alphastar=alpham[,t]-alpham[,t-1]
       alphastar[rejected]=alpham[rejected,s]
       #alphastar[rejected]=alpham[rejected,t]
       alphaused[,t]=alphastar
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
  list(rejected=rejected,decisionsm=decisionsm,cumdecisionsm=cumdecisionsm,alphaused=alphaused)
}
