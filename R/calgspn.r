# Calculate the group-sequential p-values for multiple hypotheses
calgspn=function(xm=qnorm(matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2)),
                   alpham=matrix(rep(c(0.02,0.03,0.05),each=2),ncol=3,nrow=2),
                   critm=matrix(rep(qnorm(c(0.02,0.03,0.05)),each=2),ncol=3,nrow=2),
                   matrix.list=list(diag(3),diag(3)),sided=rep(-1,2)){
  #xm: matrix of test statistics, each row is for one hypothesis and is assumed to be multivariate normal
  #alpham: matrix of alpha-levels, each row is for one hypothesis at different time points, for each row, alpha levels must be non-decreasing
  #critm: matrix of critical values, each row is for one hypothesis at different time points
  #matrix.list: list of correlation matrix corresponding to each hypothesis
  #sided: 1:  one-sided, reject if the test stat >= the critical value; 
  #       -1: one-sided, reject if the test stat <= the critical value; 
  #       0:  two-sided, reject if the absolute value of the test stat >= the critical value.
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