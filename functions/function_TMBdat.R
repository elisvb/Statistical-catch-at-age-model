TMBdat = function(dat){  
  dat.raw=readLines(dat)
  data=dat.raw[-1]
  t=sapply(data,function(x) any(grep("#",x))) 
  trues=as.vector(which(t==T))
  l=length(trues)
  list = vector("list",l)
  names(list) = gsub('\t', '', substring(names(t[t==T]),3))
  l.seq=1:l
  for(i in l.seq){
    if(i != max(l.seq)){mv=trues[i+1]-trues[i]-1} else {mv=length(data)- trues[i]}
    if(mv==1){list[[i]]=scan(dat, quiet=TRUE, skip=trues[i]+1, nlines=1)} else{
      mat = as.matrix(read.table(dat, skip=trues[i]+1, nrows=mv))
      dimnames(mat) = NULL
      list[[i]]=mat  
    } 
  }
  return(list)
}
