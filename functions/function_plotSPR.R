plotSPR=function(data=dat.,rep=rep.){
  
  RP=as.data.frame(rep$RP)
  RP$RP=c('spr20','spr30','spr40','Fmax','sprF0','F01','Flow','Fmed','Fhigh','Fmsy','SSBmsy','RECmsy','SSBmsy40','Fcol','Fmsy_noEnv','SSBmsy_noEnv','RECmsy_noEnv','Brec20','Brec40')
  names(RP)[1]='value'
  for(i in RP$RP){assign(i,RP[RP$RP==i,'value'])}
  
  pr=data.frame(f=dat$f,spr=rep$spr)
  
  yend.spr30=pr[which(abs(pr$f-spr30)==min(abs(pr$f-spr30))),'spr']
  yend.spr40=pr[which(abs(pr$f-spr40)==min(abs(pr$f-spr40))),'spr']
  
  a=ggplot(pr,aes(x=f,y=spr),environment=environment())+geom_line(size=1.5)+xlab('F')+ylab('SSB per Recruit')+
    geom_segment(aes(x=spr30, y = 0, xend = spr30, yend = yend.spr30), linetype = "dotted")+
    geom_segment(aes(x=spr40, y = 0, xend = spr40, yend = yend.spr40), linetype = "dotted")+
    geom_text(aes(x=spr30*1.05,y=yend.spr30*1.05,label=paste("F30% =",round(spr30,2))),hjust=0,vjust=1)+
    geom_text(aes(x=spr40*1.05,y=yend.spr40*1.05,label=paste("F40% =",round(spr40,2))),hjust=0,vjust=1)+
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
  
  return(a)
}
