plotYPR=function(data=dat.,rep=rep.){
  
RP=as.data.frame(rep$RP)
RP$RP=c('spr20','spr30','spr40','Fmax','sprF0','F01','Flow','Fmed','Fhigh','Fmsy','SSBmsy','RECmsy','SSBmsy40','Fcol','Fmsy_noEnv','SSBmsy_noEnv','RECmsy_noEnv','Brec20','Brec40')
names(RP)[1]='value'
for(i in RP$RP){assign(i,RP[RP$RP==i,'value'])}

pr=data.frame(f=dat$f,ypr=rep$ypr)

yend= pr[which(abs(pr$f-F01)==min(abs(pr$f-F01))),'ypr']
yend.max=pr[which(abs(pr$f-Fmax)==min(abs(pr$f-Fmax))),'ypr']
maxf=max(ypr$f)
maxypr=max(ypr$ypr)

a=ggplot(pr,aes(x=f,y=ypr))+geom_line(size=1.5)+xlab('F')+ylab('Yield per Recruit')+
  geom_segment(aes(x=F01 , y = 0, xend = F01, yend =yend), linetype = "dotted")+
  geom_text(aes(x=F01*1.1,y=yend*0.9,label=paste("F0.1 =",round(F01,2))),vjust=0,hjust=0)+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

if(Fmax>max(data$f)){
  a=a+geom_segment(aes(x=maxf*0.9,y=maxypr*0.9,xend=maxf,yend=maxypr*0.9),arrow = arrow(length = unit(0.2,"cm")))+
    geom_text(aes(x=maxf*0.8,y=maxypr*0.9),label='Fmax')
}else{
  a=a+geom_segment(aes(x=Fmax, y = 0, xend = Fmax, yend = yend.max), linetype = "dotted")+
    geom_text(aes(x=Fmax*1.1,y=yend.max*0.9,label=paste("Fmax =",round(Fmax,2))),vjust=0,hjust=0)
}
return(a)
}
