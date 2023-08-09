# TODO: adjust function arugment order to be consistent with graphical().

#' Sequential graphical procedure based on group-sequential p-values
#'
#' `seqgspgx()` implements the sequential graphical procedure for multiple
#'   hypotheses based on group-sequential p-values.
#'
#' @param pm Numeric matrix of group-sequential p-values for different
#'   hypotheses (in row) at different times (in column).
#' @param alpha Numeric scalar of the overall family-wise error rate.
#' @param W Numeric vector of the weights of the graph.
#' @param G Numeric transition matrix of the graph.
#' @return List with elements
#'   * `rejected`: the index set of rejected hypotheses
#'   * `decisionsm`: rejection decision for each endpoint (row) at each timepoint
#'     (column)
#'   * `cumdecisionsm`: cumulative rejection decision for each endpoint (row) at
#'     each timepoint (column)
#' @author Xiaodong Luo
#' @concept group-sequential p-values
#' @examples
#' seqgspgx(
#'   pm = matrix(rep(c(0.03, 0.04, 0.01), times = 2), ncol = 3, nrow = 2),
#'   alpha = 0.025,
#'   W = c(0.6, 0.4),
#'   G = rbind(c(0, 1), c(1, 0))
#' )
#' @export
seqgspgx=function(pm=matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2),alpha=0.025,
                    W=c(0.6,0.4),G=rbind(c(0,1),c(1,0))){
  n=nrow(pm)
  s=ncol(pm)
  mseq=seq(1,n,by=1)
  decisionsm=alphaused=cumdecisionsm=matrix(0,nrow=n,ncol=s)

  t=1
  abc=graphical(p=pm[,t],W=W,G=G,alpha=alpha)
  g=mseq[abc$rej.h>0]
  rejected=mseq[g]
  decisionsm[,t]=cumdecisionsm[,t]=abc$rej.h
  if (s>1){
     for (t in 2:s){
       abc=graphical(p=pm[,t],W=W,G=G,alpha=alpha)
       g=mseq[abc$rej.h>0]
       rejected=union(rejected,mseq[g])
       decisionsm[,t]=abc$rej.h
       cumdecisionsm[,t]=(abc$rej.h|cumdecisionsm[,t-1])
     }
  }
  if (length(rejected)>0) rejected=sort(rejected)
  list(rejected=rejected,decisionsm=decisionsm,cumdecisionsm=cumdecisionsm)
}
