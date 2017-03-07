plotMED=function(data=dat.,rep=rep.){
  N=exp(matrix(rep$logN,nrow=nrow(dat$propMat),byrow=FALSE))
  SR=data.frame(year=dat$y,recruit=N[1,],ssb=rep$ssb,env=dat$env1)
  SR$ssb.shift=shift.up(SR$ssb,1)

  a=ggplot(SR,aes(x=ssb.shift,y=recruit),environment=environment())+geom_point()+
    xlab('SSB')+ylab('Recruitment')+
    geom_segment(aes(x=0, y = 0, xend = max(SR$ssb)*0.95, yend = (max(SR$ssb)*0.9)/spr[which(abs(spr$f-Fmed)==min(abs(spr$f-Fmed))),'spr']), linetype = "dotted",size=1.5)+
    geom_segment(aes(x=0, y = 0, xend = max(SR$ssb)*0.95, yend = (max(SR$ssb)*0.9)/spr[which(abs(spr$f-Fhigh)==min(abs(spr$f-Fhigh))),'spr']), linetype = "dotted")+
    geom_segment(aes(x=0, y = 0, xend = max(SR$ssb)*0.95, yend = (max(SR$ssb)*0.9)/spr[which(abs(spr$f-Flow)==min(abs(spr$f-Flow))),'spr']), linetype = "dotted")+
    geom_text(aes(x=min(SR$ssb)*1.4,y=max(SR$recruit)*0.9,label=paste("Fhigh =",round(Fhigh,2))),hjust=0)+
    geom_text(aes(x=min(SR$ssb)*1.4,y=max(SR$recruit)*0.8,label=paste("Fmed =",round(Fmed,2))),hjust=0)+
    geom_text(aes(x=min(SR$ssb)*1.4,y=max(SR$recruit)*0.7,label=paste("Flow =",round(Flow,2))),hjust=0)
  print(a)
}