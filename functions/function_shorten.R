# shorten vectors/matrices/dataframes by a certain length (e.g. to make a retrospective plot)

shorten=function(x,length=1){
  if(length(dim(x))==2){
    y=x[,c(1:ncol(x)-length)]
  } else{
    y=x[c(1:length(x)-length)]
  }
  return(y)
}