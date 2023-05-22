## a function to calculate the group-sequential p-values with one endpoint
calgsp1=function(sx=qnorm(1-c(0.03,0.04,0.01)),scrit=qnorm(1-c(0.01,0.02,0.025)),salpha=c(0.01,0.02,0.025),
                     smatrix=diag(3),sided=1){
  #sx: sequence of test statistics, assumed to be multivariate normal with variance=1 and correlation matrix as smatrix
  #scrit: sequence of critical values
  #salpha: sequence of cumulative alpha levels
  #smatrix: correlation matrix of the test statistics
  #sided: 1:  one-sided, reject if the test stat >= the critical value; 
  #       -1: one-sided, reject if the test stat <= the critical value; 
  #       0:  two-sided, reject if the absolute value of the test stat >= the critical value.
  
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

