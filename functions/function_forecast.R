##########################################################################
# Function to do the forecast from TMB model
# parameters:
# x: list of the data, report and optimisation of the model (in that specific order)
# n_cast: number of time steps to predict
# n_boot: number of iterations in bootstrap
# F_cast: the way F is predicted
#             - walk: random walk for Fy such as in model
#             - zero: no fishing mortality
#             - RP: a reference point mortality (specified in RP)
#             - HCR: F changes according to Harvest Control Rule
#                     * HCR1: critical: F=0, Cautious: Fx=Fref*(SSBmsy80-SSBx)/(SSBmsy80-SSBmsy40)
# RP: value of RP mortality from model 
# data_cast: the way M (natural mortality), Weight and propMature are predicted
#             - average: constant at average last 3 years
#             - AR1age: predicted from AR1 based on age (rows)
# bpe_result: gives result of Process error forecast if True (but takes a long time!!)
# environment:
#             - average: projected as a constant value, the mean of the complete period
#             - LastYear: same as the last year of the series, allthrough the cast
#             - max: max observed throughout the series
#             - min: min observed throughout the series
#             - AR1: predicted from AR1
#             - norm: taken from normal distribution with mean and sd of the series
#
# Remarks: if I remove all rnorms, than the result is the same as the forecast in TMB, without sd
#          With rnorm, results change significantly and become much more realistic (for F=0, SSB goes up to precednings levels)
##########################################################################



forecast.mackerel=function(x=x,n_cast=10,n_boot=500,F_Cast='zero',RP='NA',data_Cast='fixed',
                           bpe_result=F,SSBmsy_result=F,environment='average'){

  require(R.utils)
  require(forecast)
  require(reshape2)
  
dat=x[[1]]
rep=x[[2]]
opt=x[[3]]
logN.ls=list()
catch.ls=list()
ssb.ls=list()
F.ls=list()
bpe.ls=list()
SSBmsy.ls=list()
ssb_casts=data.frame()
F_casts=data.frame()
rec_casts=data.frame()
catch_casts=data.frame()
SSBmsy_casts=data.frame()

# reference points
RPs=as.data.frame(rep$RP)
RPs$RP=c('spr20','spr30','spr40','Fmax','sprF0','F01','Flow','Fmed','Fhigh','Fmsy','SSBmsy','RECmsy','SSBmsy40','Fcol','Fmsy_noEnv','SSBmsy_noEnv','RECmsy_noEnv','Brec20','Brec40')
names(RPs)[1]='value'
for(s in RPs$RP){assign(s,RPs[RPs$RP==s,'value'])}

# calculate some constants (didn't report them)
phi_age = exp(opt$par['logit_ar_pe_age'])/(1 + exp(opt$par['logit_ar_pe_age']))
phi_year = exp(opt$par['logit_ar_pe_year'])/(1 + exp(opt$par['logit_ar_pe_year']))
pc_age = sqrt(1 - phi_age^2)                 
sigmape = exp(opt$par['log_std_pe'])

empty_cast=matrix(NA,ncol=n_cast,nrow=length(dat$a)) 

for(option_data_Cast in data_Cast){   
  for(option_environment in environment){
    for(option_F_Cast in F_Cast){
      if(option_F_Cast == 'RP'){RP_new=RP}else{RP_new="NoNeed"}
      for(option_RP in RP_new){
pb <- txtProgressBar(min = 0, max = 100, style = 3)        

#### Bootstrap
for(b in 1:n_boot){

## M, W, propMat (and WeightS; testing)
if(option_data_Cast=='average'){ 
  for(n in c('M','propMature','Weight','WeightS')){
    y=if(n=='WeightS') dat$ys else dat$y
    assign(paste0(n,'_cast'),cbind(dat[[n]],matrix(rep(rowMeans(dat[[n]][,((length(y)-2):length(y))]),n_cast),ncol=n_cast)))
  }
}

if(option_data_Cast=='AR1age'){
  for(n in c('M','propMature','Weight','WeightS')){
    y=if(n=='WeightS') dat$ys else dat$y
    cast=cbind(dat[[n]],empty_cast)
    for(i in 1:(length(dat$a))){
      works <- suppressWarnings(try(arima(ts(cast[1,]), order = c(1, 0, 0)), silent=TRUE))
      if ('try-error' %in% class(works) & length(unique(dat[[n]][1,]))!=1) {
        message(paste('AR1 for',n, 'failed'))
        next
      }
      if ('try-error' %in% class(works) & length(unique(dat[[n]][1,]))==1) {
        assign(paste0(n,'_cast'),cbind(dat[[n]],matrix(rep(rowMeans(dat[[n]][,((length(y)-2):length(y))]),n_cast),ncol=n_cast))) #take average
      } else {
        ari <- arima(ts(cast[1,]), order = c(1, 0, 0)) 
        fore=predict(ari,n.ahead=n_cast,se.fit = T)
        cast[i,(length(y)+1):(length(y)+n_cast)]=rnorm(n_cast,fore$pred,fore$se)
        cast[i,(length(y)+1):(length(y)+n_cast)][cast[i,(length(y)+1):(length(y)+n_cast)]<=0]=0.001
        }
    }
    assign(paste0(n,'_cast'),round(cast,digits=3))
  }
}

## environment (only env1 so far)
env_cast=c(dat$env1,rep(NA,n_cast))
if(option_environment=='min'){env_cast[(length(dat$env1)+1):length(env_cast)]=min(dat$env1)}
if(option_environment=='max'){env_cast[(length(dat$env1)+1):length(env_cast)]=max(dat$env1)}
if(option_environment=='average'){env_cast[(length(dat$env1)+1):length(env_cast)]=mean(dat$env1)}
if(option_environment=='LastYear'){env_cast[(length(dat$env1)+1):length(env_cast)]=tail(dat$env1,1)}
if(option_environment=='norm'){env_cast[(length(dat$env1)+1):length(env_cast)]=rnorm(n_cast,mean(dat$env1),sd(dat$env1))}
if(option_environment=='AR1'){ari <- arima(ts(dat$env1), order = c(1, 0, 0)) 
                              fore=predict(ari,n.ahead=n_cast,se.fit = T)
                              env_cast[(length(dat$env1)+1):length(env_cast)]=rnorm(n_cast,fore$pred,fore$se)}


## F, and zz (if no harvest rule)
F_cast=cbind(rep$'F',empty_cast)
Fy_cast=c(exp(rep$logFy),rep(NA,n_cast))
sel=exp(rep$logFa)/exp(max(rep$logFa))
for(j in 1:n_cast){
  j=j+length(dat$y)
    if(option_F_Cast=='walk') {Fy_cast[j]=exp(rnorm(1,log(Fy_cast[j-1]), exp(opt$par['logsd_logFy'])))
      if(Fy_cast[j]<0){Fy_cast[j]=0}} 
  for(i in 1:length(dat$a)){
   if(option_F_Cast=='zero') {F_cast[i,j]=0}
   if(option_F_Cast=='RP') {F_cast[i,j]=RPs[RPs$RP==option_RP,'value']*sel[i]}
  }
}
if(option_F_Cast=='walk') {F_cast=outer(exp(rep$logFa),Fy_cast)}
ZZ_cast=cbind(rep$F+dat$M,empty_cast)

# N & SSB
bpe_cast=cbind(rep$pe,empty_cast)
bpe_cast[1,]=0
logN_cast=cbind(rep$logN,empty_cast)
ssb_matrix=matrix(NA,ncol=n_cast+length(dat$y),nrow=length(dat$a))
ssb_cast=c(rep$ssb,rep(NA,n_cast))
for(j in 1:n_cast){
  j=j+length(dat$y)
  #HCR for F
  if(option_F_Cast=='HCR1') {
    if(ssb_cast[j-1]<=SSBmsy*0.4){F_cast[,j]=0}else{  # critical zone
      if(ssb_cast[j-1]>SSBmsy*0.8){F_cast[,j]=RPs[RPs$RP==option_RP,'value']*rep$Sel[,length(rep$Fy)]}else{ #Healthy zone. Fishing at indicated ref point
        F_cast[,j]=((RPs[RPs$RP==option_RP,'value']*rep$Sel[,length(rep$Fy)])*(ssb_cast[j-1]-SSBmsy*0.4))/(SSBmsy*0.8-SSBmsy*0.4)} # catious zone
    }
  }
  ZZ_cast[,j]=F_cast[,j]+M_cast[,j]
  #recruitment
  logN_cast[1,j]=rnorm(1,opt$par['R1']+log(ssb_cast[j-1])-log(1.0+exp(opt$par['R2'])*ssb_cast[j-1])+opt$par['E1']*env_cast[j-1]+0.5*exp(opt$par['logsd_logrec'])^2,exp(opt$par['logsd_logrec']))
  # PE and logN
  for(i in 2:length(dat$a)){
    mZ = phi_year*bpe_cast[i,j-1] + phi_age*(bpe_cast[i-1,j] - phi_year*bpe_cast[i-1,j-1]);
    bpe_cast[i,j]=rnorm(1,phi_year*bpe_cast[i,j-1] + phi_age*(bpe_cast[i-1,j] - phi_year*bpe_cast[i-1,j-1]),sigmape)
    logN_cast[i,j]=logN_cast[i-1,j-1]-ZZ_cast[i-1,j-1]+bpe_cast[i,j]
    if(i==length(dat$a)){
      logN_cast[i,j]=log(exp(logN_cast[i,j])+exp(logN_cast[i,j-1]+bpe_cast[i,j]-ZZ_cast[i,j-1]))
    }
  }
  ssb_cast[j]=sum(exp(logN_cast[,j])*Weight_cast[,j]*propMature_cast[,j])
  if(ssb_cast[j]<0){ssb_cast[j]=0}
}

#catches
logCpred_cast=empty_cast
cols=(length(dat$y)+1):(length(dat$y)+n_cast)
logF=F_cast
logF[logF==0]=1
logF=log(logF)
logCpred_cast[,1:n_cast]=logF[,cols]-log(ZZ_cast[,cols])+log(1-exp(-ZZ_cast[,cols]))+logN_cast[,cols]
expCpred_cast=exp(logCpred_cast)*Weight_cast[,1:n_cast]
Cpred_cast=c(rep$CpredtotW,colSums(expCpred_cast))

# msy  # does include the environment and lognormal dist!!!
if(SSBmsy_result==T){
SSBmsy_cast=c(RPs$SSBmsy,rep(NA,n_cast))
for(j in 1:n_cast){
  j=j+length(dat$y)
SSBmsy_cast[j]=SSBmsy.calc(fvector=dat$f,age.rec=min(dat$a),max.age=max(dat$a),
       wt=Weight_cast[,j],pmat=propMature_cast[,j],sel=rep$Sel[,length(dat$y)],M=M_cast[,j],bpe=bpe_cast[,j],
       alpha=exp(opt$par[which(names(opt$par)=='logR1')]),K=1/exp(opt$par[which(names(opt$par)=='logR2')]),E=opt$par[which(names(opt$par)=='E1')],
       env=env_cast[j],logsd_logrec=opt$par[which(names(opt$par)=='logsd_logrec')])
}
SSBmsy.ls[[b]]=SSBmsy_cast
}

logN.ls[[b]]=logN_cast
catch.ls[[b]]=Cpred_cast
ssb.ls[[b]]=ssb_cast
F.ls[[b]]=colMeans(F_cast)
bpe.ls[[b]]=melt(bpe_cast)

setTxtProgressBar(pb, b*100/n_boot)
}
close(pb)
cat("Averaging...")
if(bpe_result==T){bpe_casts=do.call('rbind',bpe.ls)}
if(SSBmsy_result==T){SSBmsy_casts=rbind(SSBmsy_casts,data.frame(year=seq(min(dat$y),max(dat$y)+n_cast),value=colMeans(do.call('rbind',SSBmsy.ls)),sd=apply((do.call('rbind',SSBmsy.ls)),2,sd),variable='SSBmsy',F_Cast=option_F_Cast,RP=option_RP,data_Cast=option_data_Cast,env=option_environment))}
ssb_casts=rbind(ssb_casts,data.frame(year=seq(min(dat$y),max(dat$y)+n_cast),value=colMeans(do.call('rbind',ssb.ls)),sd=apply((do.call('rbind',ssb.ls)),2,sd),variable='SSB',F_Cast=option_F_Cast,RP=option_RP,data_Cast=option_data_Cast,env=option_environment))
F_casts=rbind(F_casts,data.frame(year=seq(min(dat$y),max(dat$y)+n_cast),value=colMeans(do.call('rbind',F.ls)),sd=apply((do.call('rbind',F.ls)),2,sd),variable='F',F_Cast=option_F_Cast,RP=option_RP,data_Cast=option_data_Cast,env=option_environment))
rec_casts=rbind(rec_casts,data.frame(year=seq(min(dat$y),max(dat$y)+n_cast),value=colMeans(exp(do.call('rbind',logN.ls)[seq(1,max(dat$a)*n_boot,10),])),sd=apply(exp(do.call('rbind',logN.ls)[seq(1,max(dat$a)*n_boot,10),]),2,sd),variable='Recruitment',F_Cast=option_F_Cast,RP=option_RP,data_Cast=option_data_Cast,env=option_environment))
catch_casts=rbind(catch_casts,data.frame(year=seq(min(dat$y),max(dat$y)+n_cast),value=colMeans(do.call('rbind',catch.ls)),sd=apply((do.call('rbind',catch.ls)),2,sd),variable='Catch',F_Cast=option_F_Cast,RP=option_RP,data_Cast=option_data_Cast,env=option_environment))
      }
    }
  }
}

forecast.result=rbind(ssb_casts,F_casts,rec_casts,catch_casts,SSBmsy_casts)
if(bpe_result==T){return(list(forecast.result,bpe_casts))}else{
  return(list(forecast.result))
}
unload.Package(R.utils)
unload.Package(forecast)
unload.Package(reshape2)
}





  