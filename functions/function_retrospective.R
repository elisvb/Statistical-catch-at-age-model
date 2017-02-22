############################################################################################################################################
# Do retrospective analysis: results of same model for data set with progressively less years of data
# assumes the full model has just been run
# assumes ssb and F are in report cpp file
# outplut: list with
#           1) peel reports as dataframes
#           2) retro plot F
#           3) retro plot SSB
#           4) mohn's rho 
############################################################################################################################################

require(ggplot2)
require(plyr)

clean=function(x,ny,y,peels){
  x=as.data.frame(x)
  if(ncol(x)==ny){
    colnames(x)=dat$y
    x$age=rep(dat$a,peels)
    x$peel=rep(1:peels,each=max(dat$a))
    x=melt(x,id=c('age','peel'),variable.name='year')
    x$year=as.numeric(as.character(x$year))
    return(x)
    }else{
    colnames(x)="value"
    x$year=y
    x$peel=rep(1:peels,times=c(length(dat$y):(length(dat$y)-peels+1)))
    return(x)
  }
}

insideout <-  function(ll) {
  nms <- unique(unlist(lapply(ll, function(X) names(X))))
  ll <- lapply(ll, function(X) setNames(X[nms], nms))
  ll <- apply(do.call(rbind, ll), 2, as.list)
  lapply(ll, function(X) X[!sapply(X, is.null)])
}

retrospective=function(peels=7,file=file.,dat=dat.,para=para.,maps=maps.,rep=rep.,random=random.){
   to.cut.dat=c('y','logClower','logCupper','logCmean','crl','M','Weight','propMature','ys','S','WeightS','env1','env2')
   to.cut.para=c('logFy','logN')
   ny=length(dat$y) 
   
   ret.list=list()
   ret.list[[1]]=rep
   
   dats=dat
   y=dat$y
  for(s in 2:peels){
    dats[to.cut.dat]=lapply(dats[to.cut.dat],function(x) shorten(x))
    para[to.cut.para]=lapply(para[to.cut.para],function(x) shorten(x))
    obj <- MakeADFun(dats,para,random=random,DLL=file,map=maps)  
    opt<-nlminb(obj$par,obj$fn,obj$gr,control=list(iter.max=1000,eval.max=1000))
    ret.list[[s]]=obj$rep()
    y=c(y,dats$y)
  }
    ret=insideout(ret.list)
    ret=ret[(c('F','ssb'))]
    ret=lapply(ret, function(x) do.call('rbind.fill.matrix', x))
    ret=lapply(ret,clean,ny=ny,y=y,peels=peels)
    
    list2env(ret,envir=.GlobalEnv)

     ##plots
    ages=paste(min(dat$a),max(dat$a),sep="-")
    Fbar=ddply(F,c('year','peel'),summarise,value=mean(value))

    p1=ggplot(ssb,aes(x=year,y=value/1000,group=peel))+geom_line()+
        xlab("Year")+ylab("SSB (t)")
    p2=ggplot(Fbar,aes(x=year,y=value,group=peel))+geom_line()+
        xlab("Year")+ylab(expression(bar("F")[ages]))

    ## calculate Mohn's Rho
    ssbc=dcast(data=ssb,year~peel)
    ssbv=vector()
    Fbarc=dcast(data=Fbar,year~peel)
    Fbarv=vector()
      for(i in 2:peels){
        ro=i+1
        ssbv[i-1]=(tail(ssbc,i)[1,i]-tail(ssbc,i)[1,ro])/tail(ssbc,i)[1,i]
        Fbarv[i-1]=(tail(Fbarc,i)[1,i]-tail(Fbarc,i)[1,ro])/tail(Fbarc,i)[1,i]
      }
    mohn=data.frame(ssb=sum(ssbv),Fbar=sum(Fbarv))
    
    ## all in list
    final=list(ret,p1,p1,mohn)
    return(final)
}


