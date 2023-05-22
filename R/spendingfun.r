## Note that the OBF and Pocock spending functions are not the originally proposed ones, 
## they are the modified ones that closely resemble the original versions. 
## That being said, you might still see some differences 
spendingfun=function(alpha,fractions=seq(0.2,1,by=0.2),family="OBF",rho=1){
  #alpha: overall alpha to be spent
  #fractions: information fractions
  #family: OBF, pocock, power
  #rho: auxiliary parameter for OBF and power family
  if (family=='OBF'){
    qa=qnorm(1-alpha/2)/fractions^(rho/2)
    aseq=2*(1-pnorm(qa))
  }
  else if (family=='pocock'){
    qa=1+(exp(1)-1)*fractions
    aseq=alpha*log(qa)
  }
  else if (family=='power'){
    aseq=alpha*fractions^rho
  }
  list(aseq=aseq)
}
