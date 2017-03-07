BH.env=function(x){y=exp(unname(opt$par['R1'])+log(x)-log(1.0+exp(unname(opt$par['R2']))*x)+opt$par['E1']*SR$env.shift);return(y)}
RI.env=function(x){y=exp(unname(opt$par['R1'])+log(x)-exp(unname(opt$par['R2']))*x+opt$par['E1']*SR$env.shift);return(y)}
BH=function(x){y=exp(unname(opt$par['R1'])+log(x)-log(1.0+exp(unname(opt$par['R2']))*x));return(y)}
RI=function(x){y=exp(unname(opt$par['R1'])+log(x)-exp(unname(opt$par['R2']))*x);return(y)}


plotSR=function(data=dat.,rep=rep.,opt=opt.){
  
  RP=as.data.frame(rep$RP)
  RP$RP=c('spr20','spr30','spr40','Fmax','sprF0','F01','Flow','Fmed','Fhigh','Fmsy','SSBmsy','RECmsy','SSBmsy40','Fcol','Fmsy_noEnv','SSBmsy_noEnv','RECmsy_noEnv','Brec20','Brec40')
  names(RP)[1]='value'
  for(i in RP$RP){assign(i,RP[RP$RP==i,'value'])}
  
  N=exp(matrix(rep$logN,nrow=nrow(dat$propMat),byrow=FALSE))
  SR=data.frame(year=dat$y,recruit=N[1,],ssb=rep$ssb,env=dat$env1)
  SR$ssb.shift=shift.up(SR$ssb,1)
  SR$env.shift=shift.up(SR$env,1)
  SR$recruit.pred=exp(unname(opt$par['R1'])+log(SR$ssb.shift)-log(1.0+exp(unname(opt$par['R2']))*SR$ssb.shift)+opt$par['E1']*SR$env.shift)
  
  fnSR=if(dat$rec==5){BH}else{if(dat$rec==4){RI}else{NULL}} 
  fnSR.env=if(dat$rec==5){BH.env}else{if(dat$rec==4){RI.env}else{NULL}} 
  SR=SR[-1,]
  a=ggplot(SR,aes(x=ssb.shift,y=recruit),environment=environment())+geom_point()+
    xlab('SSB')+ylab('Recruitment')+
    stat_function(fun = fnSR.env, colour = "darkgrey")+
    stat_function(fun = fnSR, colour = "red")+
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
    geom_segment(aes(x=SSBmsy, y = 0, xend = SSBmsy, yend = BH(SSBmsy)), linetype = "dotted")+
    geom_text(x=SSBmsy*1.05, y = min(SR$recruit),label=paste0('SSBmsy = ',round(SSBmsy/1000,0),"mt"),hjust=0,vjust=0)+
    geom_segment(aes(x=0, y = 0, xend = rep$spr[dat$f==Fmsy]*max(SR$recruit)*0.20, yend = max(SR$recruit)*0.20), linetype = "dotted")+
    geom_text(x=rep$spr[dat$f==Fmsy]*max(SR$recruit)*0.20, y = max(SR$recruit)*0.20,label=paste('Fmsy =',round(Fmsy,2)),hjust=-0.1,vjust=0)+
    geom_text(x=0, y = max(SR$recruit)*0.9,label=paste('Fcol =',round(Fcol,2)),hjust=-0.1,vjust=0)+
    geom_segment(aes(x=0, y = 0, xend = max(SR$ssb.shift)*0.7, yend = max(SR$ssb.shift)*0.7/sprF0), linetype = "dotted")+
    geom_text(x=max(SR$ssb.shift)*0.7, y = max(SR$ssb.shift)*0.7/sprF0,label=paste('sprF0 =',round(sprF0,2)),hjust=-0.1,vjust=1)
  
  return(a)
}
