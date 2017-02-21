##### crlTransform #####
## continuous ratio log transform 
##   Cadigan, N. 2016a. A state-space stock assessment model for northern cod, 
##   including under-reported catches and variable natural mortality rates. 
##   Can. J. Fish. Aquat. Sci., 73: 296–308.

crlTransform=function(x){
  catch.num=t(x)
  na=ncol(catch.num) 
  ny=nrow(catch.num) 
  ctot = apply(catch.num,1,sum) 
  ctot.matrix = matrix(ctot,nrow=ny,ncol=na,byrow=FALSE) 
  catnump = catch.num
  ind = catch.num==0 
  catnump[ind]=0.5   
  ctotp = apply(catnump,1,sum) 
  ctotp.matrix = matrix(ctotp,nrow=ny,ncol=na,byrow=FALSE) 
  
  catnump[!ind] = catnump[!ind]*ctot.matrix[!ind]/ctotp.matrix[!ind] 
  pnum = catnump/ctot.matrix 
  catch.prop = t(pnum)
  cum_catch.prop = apply(catch.prop,2,cumsum)
  cond_catch.prop = catch.prop[1:(na-1),]  
  for(i in 2:(na-1)){cond_catch.prop[i,] = catch.prop[i,]/(1 - cum_catch.prop[i-1,])} 
  crl = log(cond_catch.prop/(1-cond_catch.prop))
  return(crl)
}
