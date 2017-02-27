########################################################################################
# functions for forecasting
# ASSUMES ENVIRONMENTAL EFFECT AND LOGNORMAL SR!
########################################################################################


### SPR and YPR function
spr.ypr.calc<-function(pr=pr,f,age.rec,max.age,wt,pmat,sel,M,bpe)
{
  v<-c() 
  v[1]<-1   # ratio of abundance per age over recruitment, so for first age automatically zero
  yield.a<-c()
  for(i in 2:(max.age-age.rec+1))             
  {    
    v[i]<-v[i-1]*exp(-(M[i-1]+sel[i-1]*f+bpe[i-1])) 
         if(i==max.age){v[i]=v[i]*(1+exp(-(M[i]+sel[i]*f+bpe[i])))}
  }
  SPR<-sum(v*wt*pmat)
  for(i in 1:length(v)) 
  {
    yield.a[i]<- v[i]*wt[i]*(1-exp(-sel[i]*f)) 
  }
  YPR<-sum(yield.a)
  if(pr=="spr") {return(SPR)}
  if(pr=="ypr") {return(YPR)}
} 

### SSBmsy function
SSBmsy.calc=function(fvector=fvector,age.rec=age.rec,max.age=max.age,wt=wt,pmat=pmat,sel=sel,M=M,bpe=bpe,alpha=alpha,K=K,E=E, env=env, logsd_logrec=logsd_logrec){
spr<-c()
ypr<-c()
for(i in 1:length(fvector))     
{
  spr[i]<-spr.ypr.calc(pr='spr',f=fvector[i],age.rec=age.rec,max.age=max.age,wt=wt,pmat=pmat,sel=sel,M=M,bpe=bpe) 
  ypr[i]<-spr.ypr.calc(pr='ypr',f=fvector[i],age.rec=age.rec,max.age=max.age,wt=wt,pmat=pmat,sel=sel,M=M,bpe=bpe) 
}

equi.S<-c()
equi.R<-c()
yield<-c()
for(f in 1:length(fvector))
{
  equi.S[f]<-(alpha*spr[f]*exp(E*env-(0.5*exp(logsd_logrec)*exp(logsd_logrec)))-1)*K 
  equi.R[f]<-(alpha*equi.S[f]*exp(E*env-(0.5*exp(logsd_logrec)*exp(logsd_logrec))))/(1+(equi.S[f]/K))
  yield[f]<-equi.R[f]*ypr[f]      
  if(yield[f]<0){yield[f]<-0}
} 
Fmsy<-fvector[yield==max(yield)]
Umsy<-1-exp(-Fmsy)
SSBmsy<-equi.S[yield==max(yield)]

return(SSBmsy)
}
