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

GetPrecipSites <- function(events,sites){
  i = 0
  Ns <- length(sites)
  precip.sites <- vector(mode="character",length=Nevents*Ns)
  for(season in events){
    for(site in sites){
      i=i+1
      precip.sites[i] <- paste(rundir,"precip/",season,site,"pre.in",sep="")
    }
  }
  return(precip.sites)
}

GetPrecipFiles <- function(events,sites){
  i = 0
  Ns <- length(sites)
  precip.files <- vector(mode="character",length=Nevents*Ns)
  for(season in events){
    for(site in sites){
      i=i+1
      precip.files[i] <- paste(season,site,"pre.in",sep="")
    }
  }
  return(precip.files)
}

KinerosFlowFile <- function(rundir,kin.flowfilenames,case.names,all.parameters){
  for(i in 1:Ncases){
    kinerosConn <- file(paste(rundir,kin.flowfilenames[i],sep=""))
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
    paste("  ID = 1,  LEN = 2,  WID = 0.75,  SL = ",all.parameters$lock.slope[i],",  MANNING = ",all.parameters$calib.manning[i],sep=""),
    " ",
    paste("  CV =  ",all.parameters$calib.cv[i],",   THICK = 500.  , SAT = ",all.parameters$lock.saturation[i],",  PR = 2",sep=""),
    " ",
    paste("  RELIEF =    ",all.parameters$calib.relief[i],",  SPACING =    ",all.parameters$calib.spacing[i],", IN =   ",all.parameters$calib.interception[i],",  CANOPY = ",all.parameters$lock.canopy[i],sep=""),
    " ",
    "  KS              G       DIST     POR              ROCK",
    paste("   ",all.parameters$calib.ks[i],"   ",all.parameters$calib.g[i],"  ",all.parameters$calib.distribution[i],"  ",all.parameters$lock.porosity[i],"  ",all.parameters$calib.rock[i],sep=""),
    " ",
    "  FRACT = 0.2, 0.6, 0.2    SPLASH = 50,  COH = 0.5",
    " ",
    "  Plot = H",
    " ",
    " END PLANE"), kinerosConn)
    close(kinerosConn)
  }
}

KinerosRunFile <- function(rundir,kin.runfilenames,kin.flowfilenames,precip.files,kin.fcnewfilenames,
                           kin.concfilenames,kin.outfilenames,description,duration,timestep){ 
  for(i in 1:Ncases){
    kinerosConn <- file(paste(rundir,kin.runfilenames[i],sep=""))
    writeLines(c(kin.flowfilenames[i],
                 precip.files[i],
                 kin.fcnewfilenames[i],
                 kin.concfilenames[i],
                 kin.outfilenames[i], 
                 case.names[i],
                 duration[i],
                 timestep[i],
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
 

KinerosBatchFile <- function(rundir,execdir,kin.batchmaster){
  replace.slash <- function(path.name) gsub("/","\\\\",path.name)
  kinerosConn <- file(kin.batchmaster,open="w")
    writeLines(c("cls",
                   "@echo on"),
                   kinerosConn)
    for(i in 1:Ncases){
      path.name <- paste("chdir ",execdir,case.names[i],sep="")
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
