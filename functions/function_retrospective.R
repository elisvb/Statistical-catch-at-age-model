############################################################################################################################################
# Do retrospective analysis: results of same model for data set with progressively less years of data
# assumes the full model has just been run
# outplut: plots, mohn's rho, all reports (ret, in file_retro.dat)
############################################################################################################################################
clean=function(x,ny){if(ncol(x)==ny){x=as.data.frame(x)
    colnames(x)=dat$y
    x$age=rep(dat$a,peels)
    x$peel=rep(1:peels,each=max(dat$a))
    x=melt(x,id=c('age','peel'),variable.name='year')
    x$year=as.numeric(as.character(x$year))
    return(x)
    }else{x=as.data.frame((x))
    colnames(x)="value"
    x$year=as.numeric(as.character(as.vector(rbind.fill.matrix(insideout(ret.list.dat)$y))))
    x$peel=rep(1:peels,times=c(length(dat$y):(length(dat$y)-peels+1)))
    return(x)
  }
}


retro=function(peels=7,file=file.,dat=dat.,para=para.,maps=maps.,rep=rep.,random=random.,info=''){
   to.cut.dat=c('y','logClower','logCupper','logCmean','crl','M','Weight','propMature','ys','S','WeightS','env1','env2')
   to.cut.para=c('logFy','logN')
   ny.=length(dat$y) 
   
   ret.list=list()
   ret.list.dat=list()
   ret.list[[1]]=rep
   ret.list.dat[[1]]=dat
   
  for(s in 2:peels){
    dat[to.cut.dat]=lapply(dat[to.cut.dat],function(x) shorten(x))
    para[to.cut.para]=lapply(para[to.cut.para],function(x) shorten(x))
    obj <- MakeADFun(dat2,para,random=random,DLL=file,map=maps)  
    opt<-nlminb(obj$par,obj$fn,obj$gr,control=list(iter.max=1000,eval.max=1000))
    rep=obj$rep()
    ret.list[[s]]=rep
    ret.list.dat[[s]]=dat
  }
    ret=insideout(ret.list)
    ret=ret[(c('F','logN','CpredtotW','ssb'))]
    ret=lapply(ret, function(x) do.call('rbind.fill.matrix', x))
    

    ret=lapply(ret,clean,ny=ny.)
    
    save(ret,file=paste0(file,info,"_retro.dat"))
    list2env(ret,envir=.GlobalEnv)

##make the plots
    Fbar=ddply(F,c('year','peel'),summarise,value=mean(value))
    Fbar2=ddply(F[F$age %in% 3:5,],c('year','peel'),summarise,value=mean(value))
    logN$value=exp(logN$value)
    
    # detailed one
p1=ggplot(ssb,aes(x=year,y=value/1000000,group=peel))+geom_line()+
  xlab("")+ylab("SSB")
p2=ggplot(CpredtotW,aes(x=year,y=value/1000,group=peel))+geom_line()+
  xlab("")+ylab("Predicted catches")
p3=ggplot(Fbar,aes(x=year,y=value,group=peel))+geom_line()+
  xlab("")+ylab("average F (ages 1-10)")
p4=ggplot(logN[logN$age==1,],aes(x=year,y=value,group=peel))+geom_line()+
  xlab("")+ylab("Recruitment")+
  scale_x_continuous(labels=seq(1970,2010,by=10))

png(filename=paste0(getwd(),'/IMG/retrospective/',file,info,"_Retrospective.png"),units="cm",res=200,width=16,height=11)
grid.arrange(p1,p2,p3,p4,ncol=2,bottom="Year")
dev.off()

p5=ggplot(Fbar2,aes(x=year,y=value,group=peel))+geom_line()+
  xlab("")+ylab("average F (ages 3-5)")
png(filename=paste0(getwd(),'/IMG/retrospective/',file,info,"_Retrospective_F35.png"),units="cm",res=200,width=16,height=11)
grid.arrange(p5,ncol=1,bottom="Year")
dev.off()

    #pretty and short version
p1=ggplot(ssb,aes(x=year,y=value/1000000,group=peel))+geom_line()+
  xlab("")+ylab("SSB")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(plot.margin = unit(c(0,0,0,0), "cm")) 
p2=ggplot(Fbar,aes(x=year,y=value,group=peel))+geom_line()+
  xlab("")+ylab(expression(bar("F")['1-10']))+
  theme(plot.margin = unit(c(0,0,0,0), "cm")) 

library(grid)
p1 <- ggplot_gtable(ggplot_build(p1))
p2 <- ggplot_gtable(ggplot_build(p2))
maxWidth = unit.pmax(p1$widths[2:3], p2$widths[2:3])
p1$widths[2:3] <- maxWidth
p2$widths[2:3] <- maxWidth

png(filename=paste0(getwd(),'/IMG/retrospective/',file,info,"_Retrospective_pretty.png"),units="cm",res=400,width=10,height=8)
grid.arrange(p1,p2,ncol=1,bottom="Year")
dev.off()

## calculate Mohn's Rho
ssbc=dcast(data=ssb,year~peel)
ssbv=vector()
Fbarc=dcast(data=Fbar,year~peel)
Fbarv=vector()
Fbarc2=dcast(data=Fbar2,year~peel)
Fbarv2=vector()
for(i in 2:peels){
ro=i+1
 ssbv[i-1]=(tail(ssbc,i)[1,i]-tail(ssbc,i)[1,ro])/tail(ssbc,i)[1,i]
 Fbarv[i-1]=(tail(Fbarc,i)[1,i]-tail(Fbarc,i)[1,ro])/tail(Fbarc,i)[1,i]
 Fbarv2[i-1]=(tail(Fbarc2,i)[1,i]-tail(Fbarc2,i)[1,ro])/tail(Fbarc2,i)[1,i]
}
mohn=data.frame(ssb=sum(ssbv),Fbar=sum(Fbarv),Fbar35=sum(Fbarv2))
mohn_mean=data.frame(ssb=mean(ssbv),Fbar=mean(Fbarv),Fbar35=mean(Fbarv2))
save(mohn,file=paste0("output_raw/",file,info,"mohn_sum.Rdata"))
save(mohn_mean,file=paste0("output_raw/",file,info,"mohn_mean.Rdata"))
return(mohn)
return(mohn_mean)
}


