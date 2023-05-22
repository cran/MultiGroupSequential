# Find the critical values in the general correlation matrix
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


