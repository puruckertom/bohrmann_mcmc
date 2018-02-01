#############################
#
#setup directory structure
Sys.info()

#luis in OA6
if(Sys.info()[4] == "D2626UTPURUCKE3"){
  rundir <- "c:/Dropbox/Dropbox/kineros2/"
  execdir <- "c:/Local_kineros/kineros2/"
}
#cts in OA7
if(Sys.info()[4] == "DC2626UTPURUCKE"){
  rundir <- "c:/Dropbox/kineros2/"
  execdir <- "c:/Local_kineros/kineros2/"
}
#Tommy cts in OA18
if(Sys.info()[4] == "DC2626UTBOHRMAN"){
  rundir <- "c:/Dropbox/kineros2/"
  execdir <- "c:/dec_kin/kineros2/"
}
#Tom basement PC
if(Sys.info()[4] == "PURUCKER-PC"){
  rundir <- "c:/dropbox/kineros2/"
  execdir <- "c:/Local_kineros/kineros2/"
}
#any mac, but cannot exec until we have a unix compile
if(.Platform$OS.type=="unix"){
  rundir <- path.expand("~/Dropbox/kineros2/")
  execdir <- ""
}

inputdir <- paste(rundir,"input/",sep="")
localdir <- execdir
outputdir <- paste(rundir,"output/",sep="")
datestamp <- format(Sys.time(), "%Y%b%d_%H%M")
#logfile <- paste(outputdir,"stwir_log_",datestamp,".txt",sep="")
#sink(logfile)
#showConnections(all=TRUE)
#closeAllConnections()

#subdirectories
keewokdir <- paste(execdir,"keewok/",sep="")
execdir.root <- execdir
execdir <- paste(execdir,"exec/",sep="")
#########################################


#######
R.Version()

print(paste("Started Kineros/STWIR simulation at",date()))
#install.packages("mnormt")
library(mnormt)
#install.packages("tmvtnorm")
library(tmvtnorm)
#install.packages("zoo")
library(zoo)
#install.packages("hydroGOF")
library(hydroGOF)
#install.packages("coda")
library(coda)

memory.limit()



#set some variables
Nsim = 50013 #
REM <- 0  ###  Random effects model?  1 for yes, 0 for regular model with four season variance parameters
outputdir <- paste(rundir,"output.",Nsim,"/",sep="")
if(file.exists(outputdir)==FALSE){dir.create(outputdir,showWarnings=FALSE)}
thin.coda=1
burnin = 4
start.amc = 5000
stop.amc = 12000
#ver <- "reg" lambda =1
#ver <- "log" lambda = 02
#ver <- "box.cox" Lambda = !0 !1
lambda <- .25

#load custom functions
source(paste(rundir,"01createKinerosInput.R",sep=""))

#load up all precipitation and observed flow data
source(paste(rundir,"02extractRainfallRunoff.R",sep=""))

#code to run the set of simulations
source(paste(rundir,"03aturntheCrank.R",sep=""))

#initialize loop and run kineros2 or stwir once
source(paste(rundir,"03bloopKineros.R",sep=""))

#Regular, four season variance parameters, metropolis hastings simulation
if(REM == 0) {source(paste(rundir,"04updatesMH_simplecovariance.R",sep=""))}

#adaptive metropolis hastings simulation
if(REM == 1) {source(paste(rundir,"04updatesMH.R",sep=""))}

#coda output for checking traces, etc.
source(paste(rundir,"05codaInspect.R",sep=""))
