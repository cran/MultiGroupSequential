
seqgspgx=function(pm=matrix(rep(c(0.03,0.04,0.01),times=2),ncol=3,nrow=2),alpha=0.025,
                    W=c(0.6,0.4),G=rbind(c(0,1),c(1,0))){
  ## sequential Graphical procedure based on group-sequential p-values
  #pm: group-sequential p-values
  #alpha: type-1 error rate
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