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
Ntraining <- Nsites#*training.proportion
#Ntesting <- Nsites-Ntraining
#flow.randomorder <- sample(1:Nsites)
training.sites <- flow.sites
#testing.sites <- flow.sites[flow.randomorder[(Ntraining+1):Nsites]]
#choose a season or seasons
events <- c("A","B","C","D") #("A","B","C","D")
Nevents <- length(events)
Ncases <- Nevents*Ntraining
#stitch casenames
case.names <- GetCaseFilenames(events,training.sites)
data.out <- cbind(case.names,flow.sites,days,weeks,fecal.source)
season <- substr(data.out[,1],0,1)
season <- as.factor(season)
data.out <- cbind(data.out,season)
dim(data.out)
#killD403 if it is in case names
if(length(which(case.names=="D403"))>0){
  data.out <- data.out[-which(case.names=="D403"),]
}  
dim(data.out)
write.csv(data.out,paste(rundir,"rainfall.flow.csv",sep=""))
colnames(data.out)
