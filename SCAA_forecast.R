##########################################################################
# Forecast from TMB model
##########################################################################

##### Get ready
## load everything
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

file="SCAA"

files=paste0(file,c('_dat','_rep','_opt'),'.Rdata')
out=lapply(files, function(x) get(load(x)))

## white ggplots
theme_new <- theme_set(theme_classic())
theme_new <- theme_update(axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"))

########################################################################
######## FORECAST: one scenario ########################################
########################################################################

forecast_result=forecast.mackerel(x=out,n_cast=3,n_boot=1000,
                                  F_Cast=c('RP'),RP='spr40',data_Cast='average',
                                  bpe_result=F,SSBmsy_result=F,environment=c('norm'))
forecast=forecast_result[[1]]

# PLOT: general 
forecast$max=forecast$value+forecast$sd
forecast$min=ifelse((forecast$value-forecast$sd)>0,forecast$value-forecast$sd,0)

png(filename="AllCasts.png",units="cm",res=200,width=20,height=16)
ggplot(forecast,aes(x=year,y=value))+geom_line(size=1.3)+
  facet_wrap(~variable,scale="free",ncol=2)+
  geom_ribbon(aes(ymin=min,ymax=max),alpha=0.2)+
  geom_vline(xintercept=max(out[[1]]$y),linetype='dotted')+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  xlab('Year')
dev.off()

# PLOT: health
years=1950:2023

health=forecast[forecast$variable %in% c('SSB','F'),]
health=dcast(health,formula=year~variable,value=value)
names(health)[3]='Fi'

health=health[health$year %in% years,]

ymax=max(health$Fi)*1.1
xmax=max(health$SSB)*1.1
xmin=min(health$SSB)
zones=data.frame(zone=rep(c("Critical","Cautionairy","Healthy"),each=4),x=c(xmin,out[[2]]$RP[11]*0.4,out[[2]]$RP[11]*0.4, xmin,
                                                                            out[[2]]$RP[11]*0.4,out[[2]]$RP[11]*0.8,out[[2]]$RP[11]*0.8,out[[2]]$RP[11]*0.4,
                                                                            out[[2]]$RP[11]*0.8,xmax,xmax,out[[2]]$RP[11]*0.8),
                                                                        y=c(0,0,ymax,ymax,0,0,ymax,ymax,0,0,ymax,ymax))

a=strsplit(as.character(health$year),"") 
health$year=as.numeric(unlist(ldply(a,function(x) paste(x[[3]],x[[4]],sep="")))) 

png(filename=paste0(getwd(),'/IMG/model_forecast/HarvestControl_F40_allyears.png'),units="cm",res=200,width=20,height=14)
ggplot(health,aes(x=SSB,y=Fi))+geom_text(aes(label=year),vjust=0,hjust=0,size=3)+
  geom_polygon(data=zones,aes(x=x,y=y,fill=zone),alpha=0.3)+
  geom_path(linetype='dotted')+  
  geom_vline(xintercept=out[[2]]$RP[11]*0.4,linetype='dashed')+
  geom_vline(xintercept=out[[2]]$RP[11]*0.8,linetype='dashed')+
  geom_text(x=out[[2]]$RP[11]*0.4,y=max(health$Fi)*1.1,hjust=-0.1,vjust=1.1,label='SSBmsy40%')+
  geom_text(x=out[[2]]$RP[11]*0.8,y=max(health$Fi)*1.1,hjust=-0.1,vjust=1.1,label='SSBmsy80%')+
  geom_text(x=mean(c(xmin, out[[2]]$RP[11]*0.4)),y=mean(c(ymax,0)),label='Critical',size=4)+
  geom_text(x=mean(c(out[[2]]$RP[11]*0.8, out[[2]]$RP[11]*0.4)),y=mean(c(ymax,0)),label='Cautious',size=4)+
  geom_text(x=mean(c(xmax, out[[2]]$RP[11]*0.8)),y=mean(c(ymax,0)),label='Healthy',size=4)+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  scale_fill_manual(values=c('orange','red','green'))+
  theme(legend.position="none")
dev.off()


########################################################################
######## FORECAST: multiple scenarios at once ##########################
########################################################################

fore.compare=forecast.mackerel(x=out,n_cast=3,n_boot=1000,
                  F_Cast=c('walk','zero','RP'),RP=c('F01','spr40','Fmsy'),data_Cast='average',
                  bpe_result=F,SSBmsy_result=F,environment=envi)

fore.compare=fore.compare[[1]]
fore.compare$scena=paste(fore.compare$F_Cast,fore.compare$RP,sep="_") ## paste together the options for which different scenarios used

### PLOT: general
png(filename=paste0(getwd(),'/IMG/model_forecast/',file,'.',info,'.',Fy_type,n_cast,".AllCasts_compare.png"),units="cm",res=200,width=20,height=16)
ggplot(fore.compare,aes(x=year,y=value,col=scena))+geom_line(size=1.3)+
  facet_wrap(~variable,scale='free')+
  geom_ribbon(aes(ymin=(value-sd),ymax=(value+sd)),alpha=0.2,fill='black')+
  geom_vline(xintercept=max(dat$y),linetype='dotted')+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  xlab('Year')
dev.off()


