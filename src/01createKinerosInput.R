KinerosTrainingFlowFileInput <- function(){
  for(i in 1:Ntrainingcases){
    kinerosConn <- file(paste(s.dirs.exec[i],kin.flowfilenames[i],sep=""))
    writeLines(c(" BEGIN GLOBAL",
    " ",
    "  CLEN =2,  UNITS = METRIC",
    " ",
    "  DIAMS = .005, .05, .25 ! mm",
    " ",
    "  DENSITY =2.65, 2.60, 2.60 ! g/cc",
    " ",
    "  TEMP = 15 ! deg C",
    " ",
    "  Nele = 1",
    " ",
    " END GLOBAL",
    "!------------------------------------------------",
    " BEGIN PLANE",
    " ",
    paste("  ID = 1,  LEN = 2,  WID = 0.75,  SL = ",training.sites.parameters$lock.slope[i],",  MANNING = ",training.sites.parameters$calib.manning[i],sep=""),
    " ",
    paste("  CV =  ",training.sites.parameters$calib.cv[i],",   THICK = 500.  , SAT = ",training.sites.parameters$lock.saturation[i],",  PR = 2",sep=""),
    " ",
    paste("  RELIEF =    ",training.sites.parameters$calib.relief[i],",  SPACING =    ",training.sites.parameters$calib.spacing[i],", IN =   ",
          training.sites.parameters$calib.interception[i],",  CANOPY = ",training.sites.parameters$lock.canopy[i],sep=""),
    " ",
    "  KS              G       DIST     POR              ROCK",
    paste("   ",training.sites.parameters$calib.ks[i],"   ",training.sites.parameters$calib.g[i],"  ",training.sites.parameters$calib.distribution[i],"  ",
          training.sites.parameters$lock.porosity[i],"  ",training.sites.parameters$calib.rock[i],sep=""),
    " ",
    "  FRACT = 0.2, 0.6, 0.2    SPLASH = 50,  COH = 0.5",
    " ",
    "  Plot = H",
    " ",
    " END PLANE"), kinerosConn)
    close(kinerosConn)
  }
}


GetCaseFilenames <- function(events,sites){
  i = 0
  Ns <- length(sites)
  case.names <- vector(mode="character",length=Nevents*Ns)
  for(season in events){
    for(site in sites){
      i=i+1
      case.names[i] <- paste(season,site,sep="")
    }
  }
  return(case.names)  
}

GetPrecipSites <- function(sites){
  i = 0
  Ns <- length(sites)
  precip.sites <- vector(mode="character",length=Ns)
  for(site in sites){
    i=i+1
    precip.sites[i] <- paste(localdir,"precip/",site,"pre.in",sep="")
  }
  return(precip.sites)
}

GetPrecipFiles <- function(sites){
  i = 0
  Ns <- length(sites)
  precip.files <- vector(mode="character",length=Ns)
  for(site in sites){
    i=i+1
    precip.files[i] <- paste(site,"pre.in",sep="")
  }
  return(precip.files)
}


StwirTrainingRunFile <- function(){ 
  for(i in 1:Ntrainingcases){
    kinerosConn <- file(paste(s.dirs.exec[i],kin.runfilenames[i],sep=""))
    writeLines(c(kin.flowfilenames[i],
                 training.precip.files[i],
                 kin.fcnewfilenames[i],
                 kin.concfilenames[i],
                 kin.outfilenames[i], 
                 training.case.names[i],
                 duration.training[i],
                 timestep.training[i],
                 "y","n","m","n","n"),
                 kinerosConn)
    close(kinerosConn)
  }
}

KinerosRunFileOld <- function(rundir,kin.runfilenames,kin.flowfilenames,precip.files,kin.outfilenames,
                           description,duration,timestep){ 
  for(i in 1:Ncases){
    kinerosConn <- file(paste(rundir,kin.runfilenames[i],sep=""))
    writeLines(c("n",
                 kin.flowfilenames[i],
                 precip.files[i],
                 kin.outfilenames[i], 
                 case.names[i],
                 duration[i],
                 timestep[i],
                 "n","y","n","y","n"),
                 kinerosConn)
    close(kinerosConn)
  }
}
    #repeat previous run?: n
    #Parameter file: ex1.par
    #Rainfall file: ex1.pre
    #Output file:ex1.out
    #Description:Example Run
    #Duration (min):200
    #Time step (min):1
    #Courant adjustment?: n
    #Sediment?: y
    #Multipliers?: n
    #Tabular summary?: y
    #API initializing?: n
 

KinerosTrainingBatchFile <- function(){
  replace.slash <- function(path.name) gsub("/","\\\\",path.name)
  kinerosConn <- file(kin.batchmaster,open="w")
    writeLines(c("cls",
                   "@echo on"),
                   kinerosConn)
    for(i in 1:Ntrainingcases){
      path.name <- paste("chdir ",s.dirs.exec[i],sep="")
      writeLines(c(replace.slash(path.name),
                   "STWIR_executable.exe"),
                   con=kinerosConn) 
    }
    close(kinerosConn)
}
  #@echo off 
  #c: 
  #chdir c:\dropbox\kineros2
  #kineros2 < kin_input_ex1.txt

print(paste("Finished 01KinerosInput.R at",date()))
