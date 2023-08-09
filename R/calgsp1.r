##' Calculate group-sequential p-values for one hypothesis
##'
##' `calgsp1()` calculates the group-sequential p-values for one hypothesis.
##'
##' @param sx Numeric vector of test statistics, assumed to be multivariate
##'   normal with variance 1 and correlation matrix given by `smatrix`.
##' @param scrit Numeric vector of sequece of critical values for the test
##'   statistics in `sx`. It should be computed beforehand. Must have the same
##'   length as `sx`.
##' @param salpha Numeric vector of cumulative alpha levels for the test
##'   statistics in `sx`. Must have the same length as `sx`.
##' @param smatrix Matrix with the correlation matrix of the test statistics `sx`.
##' @param sided Integer scalar indicating the side of the test:
##'   * `-1`: Reject if test statistic is smaller than or equal to the critical value (one-sided)
##'   * `1`: Reject if test statistic is greater or equal to the critical value (one-sided)
##'   * `0`: Reject if the absolute value of the test statistic is greater than the critical value (two-sided)
## TODO Should we return a vector instead of a list with just one element?
##' @return List containing the group-sequential p-values.
##' @author Xiaodong Luo
##' @concept group-sequential p-values
##' @examples
##' calgsp1(
##'   sx = qnorm(1 - c(0.03, 0.04, 0.01)),
##'   scrit = qnorm(1 - c(0.01, 0.02, 0.025)),
##'   salpha = c(0.01, 0.02, 0.025),
##'   smatrix = diag(3),
##'   sided = 1
##' )
##' @export
calgsp1=function(sx=qnorm(1-c(0.03,0.04,0.01)),scrit=qnorm(1-c(0.01,0.02,0.025)),salpha=c(0.01,0.02,0.025),
                     smatrix=diag(3),sided=1){
  nx=length(sx)
  
  aa=pa=rep(0,nx)
  if (sided==-1){x=sx;crit=scrit}
  else if (sided==1){x=-sx;crit=-scrit}
  else if (sided==0){x=abs(sx);crit=scrit}
  
  if (sided==-1|sided==1){
    aa[1]=pnorm(x[1])
    pa[1]=(x[1]<=crit[1])*aa[1]+(x[1]>crit[1])*1
    if (nx>=3){
      for (i in 2:(nx-1)){
        aa[i]=OpenMx::omxMnor(covariance=smatrix[1:i,1:i], means=rep(0,i), lbound=c(crit[1:(i-1)],-Inf), ubound=c(rep(Inf,i-1),x[i]))+salpha[i-1]
        #aa[i]=mvtnorm::pmvnorm(lower=c(crit[1:(i-1)],-Inf),upper=c(rep(Inf,i-1),x[i]),sigma=smatrix[1:i,1:i])+salpha[i-1]
        pa[i]=pa[i-1]-prod(x[1:(i-1)]>crit[1:(i-1)])*(x[i]<=crit[i])*(1-aa[i])
      }
    }
    aa[nx]=OpenMx::omxMnor(covariance=smatrix[1:nx,1:nx], means=rep(0,nx), lbound=c(crit[1:(nx-1)],-Inf), ubound=c(rep(Inf,nx-1),x[nx]))+salpha[nx-1]
    #aa[nx]=mvtnorm::pmvnorm(lower=c(crit[1:(nx-1)],-Inf),upper=c(rep(Inf,nx-1),x[nx]),sigma=smatrix[1:nx,1:nx])+salpha[nx-1]
    pa[nx]=pa[nx-1]-prod(x[1:(nx-1)]>crit[1:(nx-1)])*(1-aa[nx])
  }
  else if (sided==0){
    aa[1]=2*pnorm(-x[1])
    pa[1]=(x[1]>=crit[1])*aa[1]+(x[1]<crit[1])*1
    if (nx>=3){
      for (i in 2:(nx-1)){
        aa[i]=OpenMx::omxMnor(covariance=smatrix[1:i,1:i], means=rep(0,i), lbound=c(-crit[1:(i-1)],x[i]), ubound=c(crit[1:(i-1)],Inf))*2+salpha[i-1]
        #aa[i]=mvtnorm::pmvnorm(lower=c(-crit[1:(i-1)],x[i]),upper=c(crit[1:(i-1)],Inf),sigma=smatrix[1:i,1:i])*2+salpha[i-1]
        pa[i]=pa[i-1]-prod(x[1:(i-1)]<crit[1:(i-1)])*(x[i]>=crit[i])*(1-aa[i])
      }
    }
    aa[nx]=OpenMx::omxMnor(covariance=smatrix[1:nx,1:nx], means=rep(0,nx), lbound=c(-crit[1:(nx-1)],x[nx]), ubound=c(crit[1:(nx-1)],Inf))*2+salpha[nx-1]
    #aa[nx]=mvtnorm::pmvnorm(lower=c(-crit[1:(nx-1)],x[nx]),upper=c(crit[1:(nx-1)],Inf),sigma=smatrix[1:nx,1:nx])*2+salpha[nx-1]
    pa[nx]=pa[nx-1]-prod(x[1:(nx-1)]<crit[1:(nx-1)])*(1-aa[nx])
  }
  list(pa=pa)
}

