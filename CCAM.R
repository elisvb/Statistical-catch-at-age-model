wd="C:/Users/.../GitHub/"
setwd(wd)
library(TMB)

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "function")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir(wd)

##### COMPILE MODEL #####
file="CCAM" 
compile(paste0(file,'.cpp')) 
dyn.load(dynlib(paste0(file)))

## DATA
dat=TMBdat(dat="CCAM_data.dat")
dat$crl=crlTransform(dat$crl)

## INITIAL VALUES
para=list(
  logN=matrix(10,ncol=46,nrow=10),
  logFa=c(-2, 0, 0, 0,1,1,1,1,1,1), 
  logFy=rep(-2.5,length(dat$y)),
  logsd_logS=-0.69,
  logsd_logFy=-1.48,
  logsd_logrec=-0.24,
  logsd_crl=rep(-0.5,9),
  logrecruitMean=12,
  logQ=log(1),
  R1=1, 
  R2=-12, 
  E1=0,
  E2=0,
  log_std_pe=-1.10, 
  logit_ar_pe_age=0.2,
  logit_ar_pe_year=0.4
)

## MAPPED variables 
maps <- list( logFa=factor(c(1:8,rep(8, 2))),
              logsd_crl=factor(c(3,1,2,2,2,2,2,1,1)))
if(dat$rec %in% c(1,3)){maps$R1=factor(NA)
maps$R2=factor(NA)} 
if(dat$rec>1){maps$logrecruitMean=factor(NA)}
if(all(dat$env1==0)){maps$E1=factor(NA)}
if(all(dat$env2==0)){maps$E2=factor(NA)}

## RANDOM variables
random=c("logFy","logN")

###### STATE 1: uncensored -----
obj <- MakeADFun(dat,para,random=random,map=maps,DLL=file)  

lower <- obj$par*0-Inf
upper <- obj$par*0+Inf
lower['logQ']=log(0.1);upper['logQ']=log(2)
lower[which(names(lower) =='logFa')]=log(0.01);upper[which(names(lower) =='logFa')]=log(5)

opt1<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,control=list(iter.max=1000,eval.max=1000))
opt1

rep=obj$rep()

## STATE 2: censored-----
para.new=updateParameters(opt=opt1,rep=rep,para.old=para,random=random,map=maps)
dat[c("state")]=2

obj <- MakeADFun(dat,para.new,random=random,DLL=file,map=maps)  

opt <-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,control=list(iter.max=1000,eval.max=1000))
opt

rep=obj$rep()

sdrep<-sdreport(obj,bias.correct=T)
s=summary(sdrep, "fixed", p.value = TRUE);s=round(s,digits=2)
s

plotCheck(dat,rep)    
    
save(sdrep,file=paste0(file,"_sdrep.Rdata"))
save(rep,file=paste0(file,"_rep.Rdata"))
save(opt,file=paste0(file,"_opt.Rdata"))
save(dat,file=paste0(file,"_dat.Rdata"))
