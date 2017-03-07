plotKOBE=function(data=dat.,rep=rep.){
  
  RP=as.data.frame(rep$RP)
  RP$RP=c('spr20','spr30','spr40','Fmax','sprF0','F01','Flow','Fmed','Fhigh','Fmsy','SSBmsy','RECmsy','SSBmsy40','Fcol','Fmsy_noEnv','SSBmsy_noEnv','RECmsy_noEnv','Brec20','Brec40')
  names(RP)[1]='value'
  for(i in RP$RP){assign(i,RP[RP$RP==i,'value'])}
  
  status=data.frame(year=dat$y,Harvest_rate=colMeans(rep$F)/Fmsy,Biomass=rep$ssb/SSBmsy)
  limit=5
  status[status$Harvest_rate>limit,"Harvest_rate"]=limit
  status[status$Biomass>limit,"Biomass"]=limit
  status[status$Biomass<0,"Biomass"]=0
  zones=data.frame(zone=rep(c("Recovery","Overfishing","Fishery reduction","Lightly exploited"),each=4),x=c(0,1,1,0,0,1,1,0,1,limit,limit,1,1,limit,limit,1),y=c(0,0,1,1,1,1,limit,limit,1,1,limit,limit,0,0,1,1))
  
  a=ggplot(status,aes(x=Biomass,y=Harvest_rate),environment=environment())+
    geom_polygon(data=zones,aes(x=x,y=y,fill=zone))+
    geom_text(aes(label=year),size=2)+
    geom_path(linetype='dotted')+
    geom_hline(yintercept=1)+geom_vline(xintercept = 1)+
    ylab('F/Fmsy')+xlab('SSB/SSBmsy')+
    scale_fill_manual(values=c("lightyellow","lightgreen","orange","lightyellow"))+
    scale_y_continuous(breaks=1:limit,labels=c(1:(limit-1),paste0('>',limit)))+
    scale_x_continuous(breaks=1:limit,labels=c(1:(limit-1),paste0('>',limit)))
  
  print(a)
}