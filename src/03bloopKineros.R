#############################
########
# we are keeping site.factors, obs.flow.mat, rainfall.mat 
dim(site.factors)
dim(obs.flow.mat)
dim(rainfall.mat)
# these are all 143 by xx
# randomize sites for training versus testing
# choose a site
set.seed(19)
#split training and testing sites
Nsites <- dim(obs.flow.mat)[1]
Ntraining <- length(training.seq)
Ntesting <- Nsites-Ntraining
flow.randomorder <- sample(1:Nsites)
#training.seq <- flow.randomorder[1:Ntraining]
#testing.seq <- flow.randomorder[(Ntraining+1):Nsites]
training.sites.flow <- obs.flow.mat[training.seq,]
dim(training.sites.flow)
training.sites.precip <- rainfall.mat[training.seq,]
training.sites.factors <- site.factors[training.seq,]
#testing.seq is wrong 12/23/11, has a 144?
#testing.sites.flow <- obs.flow.mat[testing.seq,]
#dim(testing.sites.flow)
#testing.sites.precip <- rainfall.mat[testing.seq,]
#testing.sites.factors <- site.factors[testing.seq,]
#Keewok was running 0.5 except for last one
#but we will switch to 1 to speed things up a bit 
timestep <- ifelse(Nendtime<250,0.5,1)

#read in parameters from Keewok
par.in <- read.table(paste(keewokdir,"parameter_231111.in",sep=""),sep="\t")
dim(par.in)
par.in <- par.in[,1:15]
par.vec<- c("slope","manning","cv","saturation","relief","spacing","interception",
            "canopy","ks","g","porosity","rock","distribution")
list.parameters <- c("lock.slope","calib.manning","calib.cv","lock.saturation","calib.relief",
                    "calib.spacing","calib.interception","lock.canopy","calib.ks","calib.g","lock.porosity",
                    "calib.rock","calib.distribution")
colnames(par.in) <- c("id","study",list.parameters)
#D403 already dropped, nned to drop A104 for some reason
if(nrow(par.in)==144){par.in <- par.in[-3,]}
#par.in[3,]
dim(par.in)
write.csv(par.in,paste(outputdir,"keewok.par.csv",sep=""))

#get initial values for relevant season/site combinations
training.sites.parameters <- par.in[training.seq,]
dim(training.sites.parameters)
testing.sites.parameters <- par.in[testing.seq,]
dim(testing.sites.parameters)

#sort out case names
if(length(case.names)==144){case.names <- case.names[-141]}
training.case.names <- case.names[training.seq]
Ntrainingcases <- length(training.case.names)
testing.case.names <- case.names[testing.seq]
Ntestingcases <- length(testing.case.names)
training.precip.files <- precip.files[training.seq]
testing.precip.files <- precip.files[testing.seq]
training.precip.sites <- precip.sites[training.seq]
testing.precip.sites <- precip.sites[testing.seq]

#####################
### conduct the first stwir/kineros2 simulation with Keeowk initial values
#create exec directories for each case - case.names
#k.dirs.exec <- paste(execdir,"k",training.case.names,"/",sep="")
s.dirs.exec <- paste(execdir,"s",training.case.names,"/",sep="")
#for(i in 1:Ncases){if(file.exists(k.dirs.exec[i])==FALSE){dir.create(k.dirs.exec[i],showWarnings=FALSE)}}
for(i in 1:Ntrainingcases){if(file.exists(s.dirs.exec[i])==FALSE){dir.create(s.dirs.exec[i],showWarnings=FALSE)}}
#Sys.chmod(paste(execdir,case.names,"/",sep="")) #mac/linux

#create flow input file names for each directory
kin.flowfilenames <- paste(training.case.names,"_flow.par",sep="")

#write precip files
#provided by Keewok- in precip.files


Ntotaltime <- end.write
#sim.flow.arr <- array(data=NA,dim=c(Ntrainingcases,Ntotaltime,Nsim))
Simnumber=1

#create run files that point to local directory files
kin.runfilenames <- paste(training.case.names,"_input.txt",sep="")
kin.outfilenames <- paste(training.case.names,"_output.txt",sep="")
kin.concfilenames <- paste(training.case.names,"_conc_output.txt",sep="")
kin.fcnewfilenames <- paste(training.case.names,"_fcnew.par",sep="")
duration.training <- duration[training.seq]
timestep.training <- timestep[training.seq]
#write batch files
kin.batchfilenames <- paste(training.case.names,"_batch.bat",sep="")
kin.batchmaster <- paste(execdir.root,"batchmaster.bat",sep="")
#kin.results <- paste(s.dirs.exec,kin.outfilenames,sep="") # changing to conc files instead
kin.results <- paste(s.dirs.exec,kin.concfilenames,sep="")
nseff.matrix.individual <- matrix(data=NA, nrow=Nsim+1, ncol=Ntrainingcases)
nseff.matrix.seasonal <- matrix(data=NA, nrow=Nsim+1, ncol=Ntrainingcases)
nseff.matrix.all <- matrix(data=NA, nrow=Nsim+1, ncol=Ntrainingcases)

#turn the crank the first time
temp.flow.arr <- TurnTheCrank(Simnumber)
#sim.flow.arr[,,Simnumber] <- temp.flow.arr

print(paste("Finished 03bloopKineros.R at",date()))
#
