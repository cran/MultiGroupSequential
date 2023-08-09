##' Calculate group-sequential p-values for multiple hypotheses
##'
##' `calgspn()` calculates the group-sequential p-values for multiple hypotheses.
##'
##' @param xm Matrix of test statistics for hypotheses (in row) and each
##'   interim (in column).
##' @param alpham Matrix of cumulative alpha levels for the test statistics `xm`.
##'   Must have the same dimensions as `xm`. For each row, alpha levels must be
##    non-decreasing.
##' @param critm Matrix of critical values for the test statistics in `xm`. It
##'   should be computed beforehand. Must have the same dimensions as `xm`.
##' @param matrix.list List of correlation matrices corresponding to each
##'   hypothesis.
##' @param sided Integer vector indicating the side of the test:
##'   * `-1`: Reject if test statistic is smaller than or equal to the critical value (one-sided)
##'   * `1`: Reject if test statistic is greater or equal to the critical value (one-sided)
##'   * `0`: Reject if the absolute value of the test statistic is greater than the critical value (two-sided)
## TODO Should we return a vector/matrix instead of a list with just one element?
##' @return List with element `pm` containing the group-sequential p-values.
##' @author Xiaodong Luo
##' @concept group-sequential p-values
##' @examples
##' calgspn(
##'   xm = qnorm(matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2)),
##'   alpham = matrix(rep(c(0.02,0.03,0.05),each=2),ncol=3,nrow=2),
##'   critm = matrix(rep(qnorm(c(0.02,0.03,0.05)),each=2),ncol=3,nrow=2),
##'   matrix.list = list(diag(3),diag(3)),
##'   sided = rep(-1,2)
##' )
##' @export
calgspn=function(xm=qnorm(matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2)),
                   alpham=matrix(rep(c(0.02,0.03,0.05),each=2),ncol=3,nrow=2),
                   critm=matrix(rep(qnorm(c(0.02,0.03,0.05)),each=2),ncol=3,nrow=2),
                   matrix.list=list(diag(3),diag(3)),sided=rep(-1,2)){
  n=nrow(xm)
  s=ncol(xm)
  pm=xm
  
  for (j in 1:n){
    crit=critm[j,]
    dcorr=matrix.list[[j]]
    pm[j,]=calgsp1(sx=xm[j,],scrit=critm[j,],salpha=alpham[j,],smatrix=matrix.list[[j]],sided=sided[j])$pa
  }
  list(pm=pm)
}
