hochbergd=function(pvalues,alpha,epsilon=1.0e-10,precision=10){
  # Hochberg procedure that can handle different alpha's for different endpoints
  m=length(pvalues);mseq=seq(1,m,by=1)
  alpha1=pmax(alpha,epsilon)
  qvalues=pvalues/alpha1
  qvalues=round(qvalues,digits=precision)
  sq=sort(qvalues,decreasing=TRUE)
  seqc=1/(mseq)
  seqc=round(seqc,digits=precision)
  ax=(sq<=seqc)
  if (sum(ax)==0){decisions=rep(FALSE,m);sqvalue=0}
  else {istar=min(mseq[ax]);sqvalue=sq[istar];decisions=qvalues<=sqvalue}
  #list(decisions=decisions,sqvalue=sqvalue,sq=sq,seqc=seqc)
  list(decisions=decisions,sqvalue=sqvalue,sq=sq)
}