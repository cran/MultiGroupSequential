### Sequential Generalized Hochberg and Hommel procedures based on q-values
seqqvalhh=function(pm=matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2),
                 alpham=matrix(rep(c(0.02,0.03,0.05),each=2),ncol=3,nrow=2),
                 epsilon=1.0e-10,precision=10,method='Hochberg'){
  ### Sequential Generalized Hochberg and Hommel procedures based on q-values
  #pm: matrix of p-values, each row is for one hypothesis at different time points
  #alpham: matrix of alpha-levels, each row is for one hypothesis at different time points, for each row, alpha levels must be non-decreasing
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