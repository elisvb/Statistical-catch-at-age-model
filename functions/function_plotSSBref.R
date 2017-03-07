plotSSBref=function(data=dat.,rep=rep.){
  RP=as.data.frame(rep$RP)
  RP$RP=c('spr20','spr30','spr40','Fmax','sprF0','F01','Flow','Fmed','Fhigh','Fmsy','SSBmsy','RECmsy','SSBmsy40','Fcol','Fmsy_noEnv','SSBmsy_noEnv','RECmsy_noEnv','Brec20','Brec40')
  names(RP)[1]='value'
  for(i in RP$RP){assign(i,RP[RP$RP==i,'value'])}
  
  ssb=data.frame(
    year=dat$y,
    bm=rep$ssb,
    sd=sdrep$sd[which(names(sdrep$value)=="ssb")]
  )
  
  SSB.F30=spr[which(abs(spr$f-spr30)==min(abs(spr$f-spr30))),'spr']*RECmsy
  SSB.F40=spr[which(abs(spr$f-spr40)==min(abs(spr$f-spr40))),'spr']*RECmsy
  
  a=ggplot(ssb,aes(x=year,y=bm/1000))+geom_line()+geom_point(size=1)+
    geom_ribbon(aes(ymin=(bm-sd)/1000,ymax=(bm+sd)/1000),alpha=0.2)+
    scale_y_continuous(expand = c(0, 0))+scale_x_continuous(expand = c(0, 0))+
    geom_hline(yintercept=SSBmsy/1000,linetype='dotted')+geom_text(x=max(ssb$year),y=SSBmsy/1000,label=paste0('SSBmsy = ',round(SSBmsy/1000,0)),hjust=1.2,vjust=-0.1)+
    geom_hline(yintercept=Brec20/1000,linetype='dotted')+geom_text(x=max(ssb$year),y=Brec20/1000,label=paste0('SSBrec20 = ',round(Brec20/1000,0)),hjust=1.2,vjust=-0.1)+
    geom_hline(yintercept=Brec40/1000,linetype='dotted')+geom_text(x=max(ssb$year),y=Brec40/1000,label=paste0('SSBrec40 = ',round(Brec40/1000,0)),hjust=1.2,vjust=-0.1)+
    geom_hline(yintercept=SSB.F30/1000,linetype='dotted')+geom_text(x=max(ssb$year),y=SSB.F30/1000,label=paste0('SSB.F30 = ',round(SSB.F30/1000,0)),hjust=3,vjust=-0.1)+
    geom_hline(yintercept=SSB.F40/1000,linetype='dotted')+geom_text(x=max(ssb$year),y=SSB.F40/1000,label=paste0('SSB.F40 = ',round(SSB.F40/1000,0)),hjust=3,vjust=-0.1)+
    ylab("SSB")+xlab('Year')
  
  return(a)
}
