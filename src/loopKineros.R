#############################
#######
#setup directory structure
#rundir <- "c:/dropbox/kineros2/"
rundir <- path.expand("~/Dropbox/kineros2/")

#subdirectories
keewokdir <- paste(rundir,"keewok/",sep="")
execdir <- paste(rundir,"exec/",sep="")

#load custom functions
source(paste(rundir,"createKinerosInput.R",sep=""))

#############################
########
# randomize sites for training versus testing
# choose a site
set.seed(19)
training.proportion <- 1
flow.sites <- c(101,102,109,110,203,204,205,206,405,406,411,412,501,502,107,108,111,112,209,210,
                307,308,407,408,105,106,503,504,201,202,401,402,403,404,409,410)
days <- c(1,1,1,1,2,2,2,2,3,3,3,3,
          1,1,1,1,2,2,2,2,3,3,3,3,
          1,1,1,1,2,2,2,2,3,3,3,3)
days <- as.factor(days)
weeks <- c(0,0,0,0,0,0,0,0,0,0,0,0,
           1,1,1,1,1,1,1,1,1,1,1,1,
           2,2,2,2,2,2,2,2,2,2,2,2)
weeks <- as.factor(weeks)
fecal.source <- c("S","P","C","X","P","X","C","P","X","S","S","C",
                 "X","S","C","P","S","P","C","S","C","X","X","P",
                 "X","P","C","X","P","S","S","X","C","P","C","S")
fecal.source <- as.factor(fecal.source)
Nsites <- length(flow.sites)
Ntraining <- Nsites*training.proportion
Ntesting <- Nsites-Ntraining
flow.randomorder <- sample(1:Nsites)
training.sites <- flow.sites[flow.randomorder[1:Ntraining]]
testing.sites <- flow.sites[flow.randomorder[(Ntraining+1):Nsites]]
#choose a season or seasons
events <- c("A","B","C","D") #("A","B","C","D")
Nevents <- length(events)
Ncases <- Nevents*Ntraining
#stitch casenames
case.names <- GetCaseFilenames(events,training.sites)
#killD403 if it is in case names
length(case.names)
if(length(which(case.names=="D403"))>0){
  data.out <- case[-which(case.names=="D403"),]
}  

######################################
#####
# read precip data and initial values of parameter inputs
# get precip input
precip.sites <- GetPrecipSites(events,training.sites)
file.exists(precip.sites)
precip.files <- GetPrecipFiles(events,training.sites)
#scan precip files for duration of experiment
rain.info <- read.table(paste(keewokdir,"rain_info.in",sep=""),sep="\t")
colnames(rain.info) <- c("id","study","duration","total.flow")
duration <- vector(mode="numeric",length=Ncases)
case.count = 0
for(j in 1:Nevents){
  for(k in 1:Ntraining){
    case.count=case.count+1
    duration[case.count] <- rain.info[which(rain.info$id==training.sites[k]&rain.info$study==events[j]),3]
  }
}
timestep <- vector(mode="numeric",length=Ncases)
timestep <- ifelse(duration<250,0.5,1)
#read observed flow
obs.flow <- read.table(paste(keewokdir,"flow_obs.in",sep=""),sep="\t")
dim(obs.flow)
obs.flow.mat <- NULL
obs.flow.mat <- matrix(data=NA,nrow=Ncases,ncol=12)
obs.flow.vec <- vector(mode="numeric",length=12)
case.count=0
for(i in 1:Nevents){
  for(j in 1:Ntraining){
    case.count=case.count+1
    obs.flow.mat[case.count,] <- as.vector(obs.flow[which(obs.flow[,1]==training.sites[j]&obs.flow[,2]==events[i]),3:14],mode="numeric")
  }
}
rownames(obs.flow.mat) <- case.names
colnames(obs.flow.mat) <- seq(5,60,5)

#read in parameters from Keewok
par.in <- read.table(paste(keewokdir,"parameter.in",sep=""),sep="\t")
dim(par.in)
par.vec<- c("slope","manning","cv","saturation","relief","spacing","interception","canopy","ks","g","porosity","rock","distribution")
colnames(par.in) <- c("id","study",par.vec)

#get initial values for each simulation season/site combination
list.parameters <- c("lock.slope","calib.manning","calib.cv","lock.saturation","calib.relief",
                     "calib.spacing","calib.interception","lock.canopy","calib.ks","calib.g","lock.porosity",
                     "calib.rock","calib.distribution")
Nparameters <- length(list.parameters)
all.parameters <- NULL
for (i in 1:Nparameters){
  keewok <- vector(mode="numeric",length=Ncases)
  case.count = 0
  for(j in 1:Nevents){
    for(k in 1:Ntraining){
      case.count=case.count+1
      keewok[case.count] <- par.in[which(par.in$id==training.sites[k]&par.in$study==events[j]),which(colnames(par.in)==par.vec[i])]
    }
  }
  temp <- assign(list.parameters[i],keewok)
  all.parameters <- c(all.parameters,list(temp))
}
names(all.parameters) <- list.parameters
rm(list=list.parameters)
rm(temp,keewok)

#####################
### conduct the first stwir/kineros2 simulation with Keeowk initial values
#create exec directories for each case - case.names
dirs.exec <- paste(execdir,case.names,"/",sep="")
for(i in 1:Ncases){if(file.exists(dirs.exec[i])==FALSE){dir.create(dirs.exec[i])}}
#Sys.chmod(paste(execdir,case.names,"/",sep="")) #mac/linux
#write flow.par files
kin.flowfilenames <- paste(case.names,"_flow.par",sep="")
KinerosFlowFile(rundir,kin.flowfilenames,case.names,all.parameters)
#write precip files
#provided by Keewok- in precip.files
#write fcpar files
kin.fcnewfilenames <- paste(case.names,"_fcnew.par",sep="")
file.copy(paste(keewokdir,"FC_new.par",sep=""),paste(rundir,kin.fcnewfilenames,sep=""))
#create run files that point to local directory files
kin.runfilenames <- paste(case.names,"_input.txt",sep="")
kin.outfilenames <- paste(case.names,"_output.txt",sep="")
kin.concfilenames <- paste(case.names,"_conc_output.txt",sep="")
KinerosRunFile(rundir,kin.runfilenames,kin.flowfilenames,precip.files,kin.fcnewfilenames,
              kin.concfilenames,kin.outfilenames,case.names,duration,timestep)
#write batch files
kin.batchfilenames <- paste(case.names,"_batch.bat",sep="")
kin.batchmaster <- paste(rundir,"batchmaster.bat",sep="")
KinerosBatchFile(rundir,execdir,kin.batchmaster)
#distribute other input files and executable files to appropriate directories
for(i in 1:Ncases){
  file.copy(paste(rundir,kin.fcnewfilenames[i],sep=""),paste(execdir,case.names[i],"/",kin.fcnewfilenames[i],sep=""),overwrite=TRUE)
  file.copy(paste(rundir,kin.batchfilenames[i],sep=""),paste(execdir,case.names[i],"/",kin.batchfilenames[i],sep=""),overwrite=TRUE)
  file.copy(paste(rundir,kin.flowfilenames[i],sep=""),paste(execdir,case.names[i],"/",kin.flowfilenames[i],sep=""),overwrite=TRUE)
  file.copy(paste(rundir,kin.runfilenames[i],sep=""),paste(execdir,case.names[i],"/",kin.runfilenames[i],sep=""),overwrite=TRUE)
  file.copy(precip.sites[i],paste(execdir,case.names[i],"/",precip.files[i],sep=""),overwrite=TRUE)
  file.copy(paste(rundir,"/stwir/STWIR_executable.exe",sep=""),paste(execdir,case.names[i],"/STWIR_executable.exe",sep=""),overwrite=TRUE)
  file.copy(paste(keewokdir,"mult.fil",sep=""),paste(execdir,case.names[i],"/mult.fil",sep=""),overwrite=TRUE)
  file.copy(paste(execdir,case.names[i],"/",kin.runfilenames[i],sep=""),paste(execdir,case.names[i],"/","kin.fil",sep=""),overwrite=TRUE)
}

####################################
#####
#run batch files- but one case at a time to keep it easy to parallelize
system.time(
  shell(kin.batchmaster) 
)
########################
  ### collect the results for initial run
#{}

#########################
Nsim=10
### Tommy MH MCMC sampling
for(i in 1:Nsim){
  #adjust parameters as necessary
  for(j in 1:Ncases){
    #run simulation for each case
  }
  #gof
}
