##### crlInversep #####
## crl transformed data to original CAA format

crlInverse=function(x){  
  na=nrow(x)
  mnew=exp(x)/(exp(x)+1)
  mprop=mnew 
  mprop[2,] = mprop[2,]*(1-mprop[1,])
  for(i in 3:na){mprop[i,] = mnew[i,]*(1 - apply(mprop[1:(i-1),],2,sum))}
  mcum=apply(mprop,2,cumsum)
  mcum=rbind(mcum,1)
  mprop=apply(mcum,2,function(x){c(x[1],diff(x))}) 
  return(mprop)
}
