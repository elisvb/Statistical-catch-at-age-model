##### updateParameters #####
## create a new parameter object, based on model output
## assumes random effects are in reported from cpp file

`%ni%` <- Negate(`%in%`) 

updateParameters = function(opt=opt.,rep=rep.,para.old=para.old.,random=random.,map=map.){
  para.new=c(split(unname(opt$par),names(opt$par)),rep[random]) # Parameter output from report (random factors) and parameter list
  for(i in names(map)[which(names(map) %in% names(para.new))]){ # mapped not as NA (fixed para)
    para.new[[i]]=unname(opt$par[which(names(opt$par)==i)][match(as.character(map[[i]]),c(1:max(as.character(map[[i]]))))])
  }  
  for(i in names(map)[which(names(map) %ni% names(para.new))]){ # mapped as NA (unused parameters)
    para.new$i=unlist(para.old[[i]])
    names(para.new)[length(para.new)]=i
  }
  para.new=para.new[names(para.old)] # get new parameters in the right order
  if((NA %in% (names(para.old)==names(para.new)))==T){stop("A parameter is missing")}else{
    if(all(names(para.old)==names(para.new))==T){return(para.new)}else{stop("The new names are no identical to the previous")}
  }
}

