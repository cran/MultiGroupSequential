## FIXME Is `crosslist()` really needed? Why not just use base::expand.grid()?
crosslist=function(b=list(a1=c(2,3),a2=c(2,4),a3=c(0,1))){
  bn=length(names(b))
  an=1
  for (i in 1:bn){an=an*length(b[[i]])}
  amatrix=matrix(0,nrow=an,ncol=bn)
  i=0;nsa=an;nsb=1
  while (i<bn){
    i=i+1
    nsa=nsa/length(b[[i]])
    if (i>=2)nsb=nsb*length(b[[i-1]])
    amatrix[,i]=rep(rep(b[[i]],each=nsb),times=nsa)
  }
  df=data.frame(amatrix)
  colnames(df)=names(b)
  list(df=df)
}
