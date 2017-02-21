##### shorten #####
## shorten vectors and matrices by a certain length (to make a retrospective plot)
## assumes years are columns

shorten=function(x,length=1){
  if(length(dim(x))==2){
    y=x[,c(1:ncol(x)-length)]
  } else{
    y=x[c(1:length(x)-length)]
  }
  return(y)
}
