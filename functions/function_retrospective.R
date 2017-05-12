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

retro=function(peels=7,file,data,para,rep,random=NA,map=NA){
   require(ggplot2)
   require(plyr)
  
   to.cut.dat=c('y','logClower','logCupper','logCmean','crl','M','Weight','propMature','ys','S','WeightS','env1','env2')
   to.cut.para=c('logFy','logN')

   ns = c("ssb","F")
   ret = array(NA, dim=c(peels,length(ns),length(data$y)),dimnames=list(1:peels,ns,data$y))
   ret[1,'ssb',]=rep$ssb
   ret[1,'F',]=colMeans(rep$'F')

  for(s in 2:peels){
    data[to.cut.dat]=lapply(data[to.cut.dat],function(x) shorten(x))
    para[to.cut.para]=lapply(para[to.cut.para],function(x) shorten(x))
    obj <- MakeADFun(data,para,random=random,DLL=file,map=map)  
    opt<-nlminb(obj$par,obj$fn,obj$gr,control=list(iter.max=1000,eval.max=1000))
    reps=obj$rep()
    ret[s,'ssb',]=c(reps$ssb,rep(NA,s-1))
    ret[s,'F',]=c(colMeans(reps$'F'),rep(NA,s-1))
  }
   ret=as.data.frame.table(ret)
   names(ret)=c('peel','variable','year','value')

     ##plots
    ages=paste(min(data$a),max(data$a),sep="-")

    p1=ggplot(ret[ret$variable=='ssb',],aes(x=year,y=value/1000000,group=peel))+geom_line()+
        xlab("Year")+ylab("SSB ('000 t)")
    p2=ggplot(ret[ret$variable=='F',],aes(x=year,y=value,group=peel))+geom_line()+
        xlab("Year")+ylab(expression(bar("F")[ages]))

    ## calculate Mohn's Rho
    ssbc=dcast(data=ret[ret$variable=='ssb',],year~peel)
    ssbv=vector()
    Fbarc=dcast(data=ret[ret$variable=='F',],year~peel)
    Fbarv=vector()
      for(i in 2:peels){
        ro=i+1
        ssbv[i-1]=(tail(ssbc,i)[1,i]-tail(ssbc,i)[1,ro])/tail(ssbc,i)[1,i]
        Fbarv[i-1]=(tail(Fbarc,i)[1,i]-tail(Fbarc,i)[1,ro])/tail(Fbarc,i)[1,i]
      }
    mohn=data.frame(ssb=sum(ssbv),Fbar=sum(Fbarv))
    mohnmean=data.frame(ssb=mean(ssbv),Fbar=mean(Fbarv))
    
    ## all in list
    final=list(ret=ret,p1=p1,p2=p2,mohn=mohn,mohnmean=mohnmean)
    return(final)
}





