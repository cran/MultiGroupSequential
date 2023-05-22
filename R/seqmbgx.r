seqmbgx=function(xm=qnorm(matrix(rep(c(0.03,0.04,0.01),times=4),ncol=3,nrow=4)),
                         informationm=matrix(rep(c(0.4,0.8,1),each=4),ncol=3,nrow=4),
                         spending=rep("OBF",4),param.spending=rep(1,4),
                         alpha=0.025,sided=-1,
                         W=c(0.5,0.5,0,0),G=rbind(c(0,0,1,0),c(0,0,0,1),c(0,1,0,0),c(1,0,0,0)),      
                         tol=1e-10,retrospective=0){
  #xm: matrix of test statistics, assumed to be multivariate normal, each row is for one hypothesis at different time points, 
  #    for each row, alpha levels must be non-decreasing
  #informationm: information for each endpoints
  #spending: type of spending function for each endpoint
  #param.spending: parameter in the spending function
  #alpha: overall alpha
  #sided: 1:  one-sided, reject if the test stat >= the critical value; 
  #       -1: one-sided, reject if the test stat <= the critical value; 
  #       0:  two-sided, reject if the absolute value of the test stat >= the critical value.
  #W: weights of the graph
  #G: transition matrix of the graph
  #tol: tolerance level for computing the critical values
  #retrospective: 0 (default) only compares the current test statistic with the updated critical value, 
  #               1 compares all the test statistics up to the current one with the updated critical values. 
  #               Even though retrospectively looking at the values are statistically valid in terms of control the type-1 error, 
  #               not retrospectively looking at the past comparisons avoids the dilemma of retrospectively inflating the alpha level.  
  
  n=nrow(xm) #number of endpoints
  s=ncol(xm) #number of analyses (interims+final)
  mseq=seq(1,n,by=1)
  crit=matrix(0,nrow=n,ncol=s)
  Hrej=rep(0,n)
  scorr=list(n)
  
  for (n1 in 1:n){
      seq.alpha=spendingfun(alpha=alpha*W[n1],fractions=informationm[n1,],family=spending[n1],rho=param.spending[n1])$aseq
      seq.alpha=pmin(seq.alpha,rep(alpha*W[n1],s))
      scorr[[n1]]=inftocor(ir=informationm[n1,])$cor
      crit[n1,]=findcrit(salpha=seq.alpha,smatrix=scorr[[n1]],sided=sided,tol=tol)$crit.value
  }
  
  
  for (s0 in 1:s){
      aset0=aset1=mseq[Hrej>0] #to capture rejected hypotheses
      ec=1;iter=0
      while (ec==1&iter<n){
         iter=iter+1
         if (retrospective==1 & sided==-1) {ax=rowSums(as.matrix(xm[,1:s0])<=as.matrix(crit[,1:s0]))}
         else if (retrospective==0 & sided==-1) {ax=as.numeric(xm[,s0]<=crit[,s0])}
         else if (retrospective==1 & sided==1) {ax=rowSums(as.matrix(xm[,1:s0])>=as.matrix(crit[,1:s0]))}
         else if (retrospective==0 & sided==1) {ax=as.numeric(xm[,s0]>=crit[,s0])}
         else if (retrospective==1 & sided==0) {ax=rowSums(as.matrix(abs(xm[,1:s0]))>=as.matrix(crit[,1:s0]))}
         else if (retrospective==0 & sided==1) {ax=as.numeric(abs(xm[,s0])>=crit[,s0])}

        
         if (sum(ax)==0){ec=0}
         else{
            aset1=union(aset0,mseq[ax>0])
            adiff=setdiff(aset1,aset0)
            an=length(adiff)
            aset0=aset1
            if (an==0){ec=0}
            else{
                Hrej[adiff]=s0;
                gdiff=setdiff(mseq,aset0)
                ng=length(gdiff)
                if (ng==0){ec=0}
                else {
                   gdiff=sort(gdiff)
                   abc=updategraph(S1=gdiff,W0=W,G0=G,S0=mseq)
                   galpha=alpha*abc$W1
                   for (n1 in 1:ng){
                        seq.alpha=spendingfun(alpha=galpha[n1],fractions=informationm[gdiff[n1],],
                                    family=spending[gdiff[n1]],rho=param.spending[gdiff[n1]])$aseq
                        seq.alpha=pmin(seq.alpha,rep(galpha[n1],s))
                        crit[gdiff[n1],]=findcrit(salpha=seq.alpha,smatrix=scorr[[gdiff[n1]]],sided=sided,tol=tol)$crit.value
                   }
                }
            }
         }
      }
  }
  rejected=NULL
  if (sum(Hrej)>0) rejected=sort(mseq[Hrej>0])
  cumdecisionsm=matrix(0,nrow=n,ncol=s)
  for (s0 in 1:s){
    cumdecisionsm[Hrej<=s0&Hrej>=1,s0]=1
  }
  decisionsm=cumdecisionsm[,s]
  
  list(Hrej=Hrej,rejected=rejected,decisionsm=decisionsm,cumdecisionsm=cumdecisionsm)
}
