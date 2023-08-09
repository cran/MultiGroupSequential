##' Check critical values
##'
##' `checkcrit()` is a helper function that checks if the critical values
##' are valid.
##'
##' @param scrit Numeric vector of critical values.
##' @param salpha Numeric vector of cumulative alpha levels.
##' @param smatrix General correlation matrix.
##' @inherit calgspn
## TODO Add further details about critical values returned (`sc` in the function)
##' @return List with:
##'   - `crit.value`: Critical values
##'   - `salpha`: Cumulative alpha levels passed to `salpha` argument
##' @author Xiaodong Luo
## TODO Do we really need this in such small and specialized package?
##' @concept group-sequential
##' @concept critical values
##' @concept efficacy boundary
##' @examples
##' checkcrit(
##'   scrit = qnorm(c(0.01, 0.02, 0.025)),
##'   salpha = c(0.01, 0.02, 0.025),
##'   smatrix = diag(3),
##'   sided = 1
##' )
##' @export
checkcrit=function(scrit=qnorm(c(0.01,0.02,0.025)),salpha=c(0.01,0.02,0.025),smatrix=diag(3),sided=1){
  #scrit: sequence of critical values
  #salpha: sequence of cumulative alpha levels, to be compared with computed alpha levels
  #smatrix: general correlation matrix
  #sided: 1:  one-sided, reject if the test stat >= the critical value; 
  #       -1: one-sided, reject if the test stat <= the critical value; 
  #       0:  two-sided, reject if the absolute value of the test stat >= the critical value.
  
  ns=length(scrit)
  
  sc=rep(0,ns)
  if (sided==1)sc[1]=1-pnorm(scrit[1])
  else if (sided==0)sc[1]=2*(1-pnorm(scrit[1]))
  else if (sided==-1)sc[1]=pnorm(scrit[1])
  
  if (ns>=2){
    for (j in 2:ns){
      if (sided==0){
        sc[j]=OpenMx::omxMnor(covariance=smatrix[1:j,1:j], means=rep(0,j), lbound=c(-scrit[1:(j-1)],scrit[j]), ubound=c(scrit[1:(j-1)],Inf))*2+sc[j-1]
      }
      else if (sided==1){
        sc[j]=OpenMx::omxMnor(covariance=smatrix[1:j,1:j], means=rep(0,j), lbound=c(rep(-Inf,j-1),scrit[j]), ubound=c(scrit[1:(j-1)],Inf))+sc[j-1]
      }
      else if (sided==-1){
        sc[j]=OpenMx::omxMnor(covariance=smatrix[1:j,1:j], means=rep(0,j), lbound=c(scrit[1:(j-1)],-Inf), ubound=c(rep(Inf,j-1),scrit[j]))+sc[j-1]
      }
    }
  }
  list(alpha.level=sc,salpha=salpha)
}
