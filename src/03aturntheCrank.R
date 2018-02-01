
###
TurnTheCrank <- function(Simnumber){
  ########################
  #write kineros flow parameter input files to exec directories
  KinerosTrainingFlowFileInput()
  
  #write fcpar files

  file.copy(paste(keewokdir,"FC_new.par",sep=""),paste(s.dirs.exec,kin.fcnewfilenames,sep=""),overwrite=TRUE)
  file.exists(paste(s.dirs.exec,kin.fcnewfilenames,sep=""))
  
  #check this for dimension errors- should be 107 checked 11/23/11
  StwirTrainingRunFile()

  KinerosTrainingBatchFile()
  #distribute other input files and executable files to appropriate directories
  if(Simnumber==1){
    for(i in 1:Ntrainingcases){
      #file.copy(paste(execdir.root,kin.fcnewfilenames[i],sep=""),paste(s.dirs.exec[i],kin.fcnewfilenames[i],sep=""),overwrite=TRUE)
      #file.copy(paste(execdir.root,kin.batchfilenames[i],sep=""),paste(s.dirs.exec[i],kin.batchfilenames[i],sep=""),overwrite=TRUE)
      #file.copy(paste(execdir.root,kin.flowfilenames[i],sep=""),paste(s.dirs.exec[i],kin.flowfilenames[i],sep=""),overwrite=TRUE)
      #file.copy(paste(execdir.root,kin.runfilenames[i],sep=""),paste(s.dirs.exec[i],kin.runfilenames[i],sep=""),overwrite=TRUE)
      file.copy(training.precip.sites[i],paste(s.dirs.exec[i],training.precip.files[i],sep=""),overwrite=TRUE)
      file.copy(paste(execdir.root,"/stwir/STWIR_executable.exe",sep=""),paste(s.dirs.exec[i],"STWIR_executable.exe",sep=""),overwrite=TRUE)
      file.copy(paste(keewokdir,"mult.fil",sep=""),paste(s.dirs.exec[i],"mult.fil",sep=""),overwrite=TRUE)
      file.copy(paste(s.dirs.exec[i],kin.runfilenames[i],sep=""),paste(s.dirs.exec[i],"kin.fil",sep=""),overwrite=TRUE)
    }
  }
  ####################################
  #####
  #STWIR run batch files- but one case at a time to keep it easy to parallelize
  elapsed.time <- system.time(
    shell(kin.batchmaster) 
  )
  print(paste(round(elapsed.time[3]),"seconds elapsed"))
  
  
  ########################
  #get the relevant observed flow values and create an observed array with sim/training dimensions
  Ntotaltime <- end.write
  obs.flow.arr <-  array(data=NA,dim=c(Ntrainingcases,Ntotaltime))  
  dim(obs.flow.mat)
  #View(obs.flow.mat[,219:279])
  dim(obs.flow.arr)
  obs.flow.arr <- obs.flow.mat[training.seq,]
  #View(obs.flow.arr[,219:279])  
  #summary(obs.flow.arr)
  dim(obs.flow.arr)
  
  #View(obs.flow.arr[,209:279])
  #write.csv()
  ### collect the results for initial run
  sim.flow.arr <- array(data=NA,dim=c(Ntrainingcases,Ntotaltime))
  dim(sim.flow.arr) #107 279 Nsim
  #rownames and colnames to help with qa
  rownames(sim.flow.arr) <- rownames(obs.flow.arr)
  colnames(sim.flow.arr) <- colnames(obs.flow.mat)
  Nendtime.training <- Nendtime[training.seq]
  timestep.training <- timestep[training.seq]
  Nreadlines <- Nendtime.training/timestep.training+1
  length(Nreadlines) #107
  length(kin.results) #107
  #report back the simulation results
  for(i in 1:Ntrainingcases){  #Ntrainingcases
    #temp.skip <- read.table(file=kin.results[i],fill=TRUE,colClasses = "character",blank.lines.skip=FALSE,col.names=paste("c",1:12,sep=""))
    #View(temp.skip)
    #Ntoskip <- which(temp.skip[,1]=="Elapsed") + 2 #changed to conc output
    Ntoskip <- 3
    #temp.sim <- scan(file=kin.results[i],skip=Ntoskip,what=list("","","",""),nlines=Nreadlines[i])
    temp.sim <- read.table(file=kin.results[i],skip=Ntoskip,header=TRUE)
    class(temp.sim$X.min.)
    class(temp.sim$X.mm.)
    temp.mat <- rbind(temp.sim$X.min.,temp.sim$X.mm.)
    #View(temp.mat)
    class(temp.mat)
    dim(temp.mat)
    #drop any timesteps not at integer resolution
    temp.seq <- 1:Nendtime.training[i]
    print(paste(i,rownames(sim.flow.arr)[i],Nendtime.training[i]))
    print(paste("Nreadlines=", Nreadlines[i]))
    print(paste("Nendtime=", Nendtime.training[i]))
    print(paste("dim temp.mat=",dim(temp.mat)))
    print(paste("first row =",temp.mat[,1]))
    print(paste("last row =",temp.mat[,dim(temp.mat)[[2]]]))
    #extract them
    temp.sim.vec <- temp.mat[2,which(temp.mat[1,]%in%temp.seq)]
    class(temp.sim.vec)
    Nsimtime <- length(temp.sim.vec)
    print(paste("Nsimtime=", Nsimtime))
    Nstartsimtime <- Ntotaltime-Nsimtime
    #paste them into sim object
    sim.flow.arr[i,(Nstartsimtime+1):Ntotaltime] <- temp.sim.vec   
    print(paste("loop:",Simnumber,"; case:",i))
  }
  
  #kluge to fix NA problem at minute 60/279
  #sim.flow.arr[which(is.na(sim.flow.arr[,279,Simnumber])),279,Simnumber] <- sim.flow.arr[which(is.na(sim.flow.arr[,279,Simnumber])),278,Simnumber]
  #sim.flow.arr[,279,Simnumber]
  
  #convert deltas to cumulative
  #dim(sim.flow.arr)
  #dim(obs.flow.arr)
  #sim.flow.cumulative.arr <- array(data=NA,dim=dim(sim.flow.arr))
  #dim(sim.flow.cumulative.arr)
  #obs.flow.cumulative.arr <- array(data=NA,dim=dim(obs.flow.arr))
  #for(i in 1:Ntrainingcases){
  #  for(j in 1:Ntotaltime){
  #    if(!is.na(sim.flow.arr[i,j,1])){
  #      sim.flow.cumulative.arr[i,j,1] <- sum(sim.flow.arr[i,1:j,Simnumber],na.rm=TRUE)
  #    }
  #    if(!is.na(obs.flow.arr[i,j])){
  #      obs.flow.cumulative.arr[i,j] <- sum(obs.flow.arr[i,1:j],na.rm=TRUE)
  #    }
  #  }
  #}
  
  #plot the predicted versus the observed time series for each plot
  dim(sim.flow.arr)
  dim(obs.flow.arr)
  rownames(sim.flow.arr) <- rownames(obs.flow.arr)
  colnames(sim.flow.arr) <- colnames(obs.flow.mat)
  #View(sim.flow.arr[,219:279,1])
  
  #calculate seasonal and overall means for nash sutcliffe
  stopA <- ceiling(dim(obs.flow.arr)[[1]]*.25)
  stopB <- ceiling(dim(obs.flow.arr)[[1]]*.25)*2
  stopC <- ceiling(dim(obs.flow.arr)[[1]]*.25)*3
  stopD <- dim(obs.flow.arr)[[1]]
  
  mean.obs.A <- mean(obs.flow.arr[1:stopA,],na.rm=TRUE)
  mean.obs.B <- mean(obs.flow.arr[(stopA+1):stopB,],na.rm=TRUE)
  mean.obs.C <- mean(obs.flow.arr[(stopB+1):stopC,],na.rm=TRUE)
  mean.obs.D <- mean(obs.flow.arr[(stopC+1):stopD,],na.rm=TRUE)
  mean.obs.All <- mean(obs.flow.arr,na.rm=TRUE)
  
  #produce nseff values for individual, seasonal, and global means
  nseff.individual <- vector(mode = "numeric", length = Ntrainingcases)
  nseff.seasonal <- vector(mode = "numeric", length = Ntrainingcases)
  nseff.all <- vector(mode = "numeric", length = Ntrainingcases)
  for(i in 1:Ntrainingcases){
    if(Simnumber == 1 | Simnumber == 2){
      gof.sim <- as.vector(sim.flow.arr[i,])
    }else{
      gof.sim <- as.vector(accepted.flow.arr[i,])
    }
    #change to training.seq for observed
    gof.obs <- as.vector(obs.flow.arr[i,])
    tsbind <- matrix(rbind(gof.obs,gof.sim),nrow=2,ncol=279)
    nseff.individual[i] <- suppressWarnings(NSE(gof.sim,gof.obs))
    obs.indices <- as.vector(which(!is.na(obs.flow.arr[i,])))
    Nobs <- length(obs.indices)
    NSE.numerator <- sum((gof.obs[obs.indices]-gof.sim[obs.indices])^2,na.rm=TRUE)
    if(i <= stopA){
      NSE.denominator.seasonal <- sum((gof.obs[obs.indices]-mean.obs.A)^2)
    }else if(i <= stopB){
      NSE.denominator.seasonal <- sum((gof.obs[obs.indices]-mean.obs.B)^2)
    }else if(i <= stopC){
      NSE.denominator.seasonal <- sum((gof.obs[obs.indices]-mean.obs.C)^2)
    }else{
      NSE.denominator.seasonal <- sum((gof.obs[obs.indices]-mean.obs.D)^2)
    }
    NSE.denominator.all <- sum((gof.obs[obs.indices]-mean.obs.All)^2)
    nseff.seasonal[i] <- round(1-(NSE.numerator/NSE.denominator.seasonal),digits=2)
    nseff.all[i] <- round(1-(NSE.numerator/NSE.denominator.all),digits=2)
  }
  
  prop.nseff.individual <- round(length(which(nseff.individual>0))/length(nseff.individual),digits=2)
  prop.nseff.seasonal <- round(length(which(nseff.seasonal>0))/length(nseff.seasonal),digits=2)
  prop.nseff.all <- round(length(which(nseff.all>0))/length(nseff.all),digits=2)
  if(Simnumber <= Nsim){  
    nseff.matrix.individual[Simnumber,] <- nseff.individual
    nseff.matrix.seasonal[Simnumber,] <- nseff.seasonal
    nseff.matrix.all[Simnumber,] <- nseff.all
  }else{
    nseff.mean.individual <<- nseff.individual
    nseff.mean.seasonal <<- nseff.seasonal
    nseff.mean.all <<- nseff.all    
  }
  
  if(Simnumber %% 100 == 0 | Simnumber == 1 | Simnumber ==2 | Simnumber == 10 | Simnumber == Nsim+1){ 
    #produce nseff plots versus seasonal means
    pdf(paste(outputdir,Nsim,".nseff.",Simnumber,".ts.pdf",sep=""),width=10.5,height=7.5)
      par(mfrow=c(3,3))
      prop.nseff.individual.local <- round(length(which(nseff.individual>0))/length(nseff.individual),digits=3)
      prop.nseff.seasonal.local <- round(length(which(nseff.seasonal>0))/length(nseff.seasonal),digits=3)
      prop.nseff.all.local <- round(length(which(nseff.all>0))/length(nseff.all),digits=3)
      hist(nseff.individual,col="cornflowerblue",main=paste("p(nseff>0)=",prop.nseff.individual.local))
      hist(nseff.seasonal,col="cornflowerblue",main=paste("p(nseff>0)=",prop.nseff.seasonal.local))
      hist(nseff.all,col="cornflowerblue",main=paste("p(nseff>0)=",prop.nseff.all.local))    
      for(i in 1:Ntrainingcases){
        if(Simnumber == 1 | Simnumber == 2 | Simnumber == Nsim + 1){
          gof.sim <- as.vector(sim.flow.arr[i,])
        }else{
          gof.sim <- as.vector(accepted.flow.arr[i,])
        }
        #change to training.seq for observed
        gof.obs <- as.vector(obs.flow.arr[i,])
        tsbind <- matrix(rbind(gof.obs,gof.sim),nrow=2,ncol=279)
        compare.flow.zoo <- zoo(t(tsbind),order.by=sim.time)
        title.main <- paste("NSE_I=",nseff.individual[i],"NSE_S=",nseff.seasonal[i],"NSE_A=",nseff.all[i],"\n",
                            "Sim=",rownames(sim.flow.arr)[i],"Flow=",rownames(obs.flow.arr)[i]) #,"\n",last.warning
        ifelse(nseff.seasonal[i]>0,pl.col <- c(1,3),pl.col <- c(1,2))
        plot.zoo(compare.flow.zoo,type="b",pch=1:2,col=pl.col,plot.type="single",
                 ylab="Flow",xlab="Experiment Time",main=title.main)
      }
    dev.off()
  }

  #or comment that out and return cumulative
  #sim.flow.cumulative.arr[,,Simnumber]
  return(sim.flow.arr[,])
}
#wtf

