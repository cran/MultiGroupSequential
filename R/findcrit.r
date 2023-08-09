##' Calculate critical values
##'
##' `findcirt()` calculates the critical values in the general correlation matrix
##'
## TODO Add further details about how these critical values are computed
##' @inherit checkcrit
##' @param tol Numeric scalar with the tolerance level for computing
##'   critical values.
## TODO Explain `alpha.tol` further
##' @param alpha.tol Numeric scalar. If the alpha increment is less than this,
##'   the critical value is set to a large number determined by `alpha.tol`.
## TODO Why not just return `sc` instead of `list(crit.value = sc)`?
##' @return List with element `crit.value` containing the obtained critical values.
##' @author Xiaodong Luo
##' @concept group-sequential
##' @concept critical values
##' @concept efficacy boundary
##' @examples
##' findcrit(
##'   salpha = c(0.01, 0.02, 0.025),
##'   smatrix = diag(3),
##'   sided = 1,
##'   tol = 1e-10,
##'   alpha.tol = 1e-11
##' )
##' @export
## TODO Refactor multiple calls to `omxMnor()` to make it simpler?
findcrit=function(salpha=c(0.01,0.02,0.025),smatrix=diag(3),sided=1,tol=1e-10,alpha.tol=1e-11){
  #salpha: sequence of cumulative alpha levels
  #smatrix: general correlation matrix
  #sided: 1:  one-sided, reject if the test stat >= the critical value; 
  #       -1: one-sided, reject if the test stat <= the critical value; 
  #       0:  two-sided, reject if the absolute value of the test stat >= the critical value.
  #tol: tolerance level for computing the critical values
  #alpha.tol: if the alpha increment is less than this, then the critical value is set to a large number determined by alpha.tol
  
  ns=length(salpha)
  
  sc=rep(0,ns)
  if (sided==1)sc[1]=qnorm(1-salpha[1])
  else if (sided==0)sc[1]=qnorm(1-salpha[1]/2)
  else if (sided==-1)sc[1]=qnorm(salpha[1])
  
  
  if (ns>=2){
    for (j in 2:ns){
      if (salpha[j]-salpha[j-1]>alpha.tol){
        xfunc=function(x,xcrit=rep(0,j-1)){
          if (sided==0){
            temp=OpenMx::omxMnor(covariance=smatrix[1:j,1:j], means=rep(0,j), lbound=c(-xcrit,x), ubound=c(xcrit,Inf))*2-(salpha[j]-salpha[j-1])
          }
          else if (sided==1){
            temp=OpenMx::omxMnor(covariance=smatrix[1:j,1:j], means=rep(0,j), lbound=c(rep(-Inf,j-1),x), ubound=c(xcrit,Inf))-(salpha[j]-salpha[j-1])
          }
          else if (sided==-1){
            temp=OpenMx::omxMnor(covariance=smatrix[1:j,1:j], means=rep(0,j), lbound=c(xcrit,-Inf), ubound=c(rep(Inf,j-1),x))-(salpha[j]-salpha[j-1])
          }
          temp
        }
        if (sided==0) tinterval=c(0,qnorm(1-alpha.tol/2))
        else if (sided==1) tinterval=c(0,qnorm(1-alpha.tol))
        else if (sided==-1) tinterval=c(-qnorm(1-alpha.tol),0)
        sc[j]<- uniroot(xfunc, interval=tinterval, tol = tol, xcrit=sc[1:(j-1)])$root
      }
      else if (salpha[j]-salpha[j-1]<=alpha.tol){
        sc[j]=(sided==0)*qnorm(1-alpha.tol/2)+(sided==1)*qnorm(1-alpha.tol)+(sided==-1)*qnorm(alpha.tol)
      }
    }
  }
  list(crit.value=sc)
}


