#############################
########
# randomize sites for training versus testing
# choose a site

training.proportion <- .75
flow.sites <- c(101,102,109,110,203,204,205,206,405,406,411,412,501,502,107,108,111,112,209,210,
                307,308,407,408,105,106,503,504,201,202,401,402,403,404,409,410)
all.sites.unique <- sort(flow.sites)
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

#we can only sort the cbind
factors.df <- cbind(flow.sites,days,weeks,fecal.source)
head(factors.df)
factors.sorted.df <- factors.df[order(flow.sites),]
#View(factors.sorted.df)

factors.144 <- rbind(factors.sorted.df,factors.sorted.df,factors.sorted.df,factors.sorted.df)
dim(factors.144)
head(factors.144)

# wet or dry
wet.dry <- c(rep("wet",72),rep("dry",72)) 
factors.144 <- cbind(factors.144,wet.dry)
head(factors.144)

#seasons
seasons4 <- c(rep("A",36),rep("B",36),rep("C",36),rep("D",36))
factors.144 <- cbind(seasons4,factors.144)
head(factors.144)
factors.144[,1] <- as.factor(factors.144[,1])

#sort (order) ascending
class(factors.144)
factors.144[order(factors.144[,1],factors.144[,2]),]

factors.143 <- factors.144[-131,]

#select sites for training
Nsites <- 36
Ntraining <- Nsites*training.proportion #27
Ntrainingunique <- Ntraining
Ntesting <- Nsites-Ntraining #9
flow.sites <- flow.sites[sort.list(flow.sites)]
all.sites <- paste(seasons4,flow.sites,sep="")

# have not dealt with d403 here, depends on who gets it (randomly)
all.seq <- 1:36
set.seed(6)
training.seq <- sample(1:36,Ntraining,replace=FALSE)
training.seq <- sort(training.seq)
training.seq.unique <- training.seq
training.sites.unique <- flow.sites[training.seq.unique]
testing.seq <- all.seq[-training.seq]
testing.seq <- sort(testing.seq)
testing.seq.unique <- testing.seq
testing.sites.unique <- flow.sites[testing.seq.unique]
head(training.seq)

whereD403 <- NULL
wherenotD403 <- NULL

training.seq <- c(training.seq,training.seq+36,training.seq+72,training.seq+108)
training.sites <- all.sites[training.seq]
"D403" %in% training.sites
length(training.seq)
length(training.sites)
training.sites

if("D403" %in% training.sites == TRUE){
  print("true")
  length(training.seq)
  length(training.sites)
  whereD403 <- which(training.sites=="D403")
  whereD403
  training.seq <- training.seq[-whereD403]
  length(training.seq)
  training.sites <- training.sites[-whereD403]
  length(training.sites)
  training.seq #this is still wrong because it goes to 144
  training.sites#but this is right
  #update training.seq to work with 143 length inputs
  all.sites[training.seq]
  training.sites
  training.seq[whereD403:length(training.seq)]
  training.seq[whereD403:length(training.seq)] = training.seq[whereD403:length(training.seq)]-1
}else{
  print("false")
  training.sites.unique
  wherenotD403 <- 3*length(training.sites.unique)+(min(which(training.sites.unique>"403")))
  wherenotD403
  training.sites.unique[wherenotD403-3*length(training.sites.unique)]
  training.seq[1:(wherenotD403-1)]
  #
  training.seq[wherenotD403:length(training.sites)]
  training.seq[wherenotD403:length(training.sites)] <- training.seq[wherenotD403:length(training.sites)]-1
  training.seq[wherenotD403:length(training.sites)]
  training.sites
}
cbind(training.sites,training.seq)

#this is 143 long

training.seq
length(training.seq) == length(unique(training.seq))
factors.143[training.seq,]
"D403" %in% training.sites
whereD403
wherenotD403
length(training.seq)

testing.seq <- c(testing.seq,testing.seq+36,testing.seq+72,testing.seq+108)
testing.sites <- all.sites[testing.seq]
if("D403" %in% testing.sites){
  length(testing.sites)
  whereD403 <- which(testing.sites=="D403")
  testing.seq <- testing.seq[-whereD403]
  testing.sites <- all.sites[testing.seq]
  testing.sites
  testing.seq[whereD403:length(testing.seq)]
  testing.seq[whereD403:length(testing.seq)] = testing.seq[whereD403:length(testing.seq)]-1
  length(testing.sites)
}

Ncases <- length(all.sites)

#flow.randomorder <- sample(1:Nsites)
#training.sites <- flow.sites
#testing.sites <- flow.sites[flow.randomorder[(Ntraining+1):Nsites]]

#choose a season or seasons
#events <- c("A","B","C","D") #("A","B","C","D")
#Nevents <- length(events)
#Ncases <- Nevents*Ntraining

#stitch casenames
#deprecated
case.names <- all.sites

#data.out <- cbind(case.names,flow.sites,days,weeks,fecal.source,wet.dry)
data.out <- factors.144
season <- seasons4
season <- as.factor(season)
data.out <- cbind(data.out,season)
dim(data.out)

#killD403 if it is in case names
if(length(which(case.names=="D403"))>0){
  Ncases=Ncases-1
  data.out <- data.out[-which(case.names=="D403"),]
  case.names <- all.sites[-which(case.names=="D403")]
}  
dim(data.out)
#site.factors contains all relevant factors for the 144 sampling events
site.factors <- data.out
write.csv(data.out,paste(execdir.root,"rainfall.flow.csv",sep=""))
colnames(data.out)
length(case.names)

######################################
#####
# read precip data and initial values of parameter inputs
# get precip input
precip.sites <- GetPrecipSites(all.sites)
killD403 <- which(file.exists(precip.sites)==FALSE)
precip.files <- GetPrecipFiles(all.sites)
#killd403 again
if(length(precip.sites)==144){
  precip.sites <- precip.sites[-killD403]
  precip.files <- precip.files[-killD403]
  #case.names <- case.names[-131]
}

length(precip.sites)
Nprecip <- length(precip.sites)
Ntimesteps <- vector(mode="numeric",length=Nprecip)
Nmaxprecip <- vector(mode="numeric",length=Nprecip)
Ntimesteps <- vector(mode="numeric",length=Nprecip)
Nendtime <- vector(mode="numeric",length=Nprecip)

#find max length of time series and other associated values
for(i in 1:Nprecip){
  precip.data <- scan(precip.sites[i],skip=7,what=c("numeric","numeric"),na.strings="END")
  precip.data <- as.vector(as.numeric(precip.data[-length(precip.data)]))
  dim(precip.data) <- c(2,length(precip.data)/2)    
  Ntimesteps[i] <- dim(precip.data)[2]
  Nmaxprecip[i] <- precip.data[2,(Ntimesteps[i]-1)]
  Nendtime[i] <- precip.data[1,(Ntimesteps[i]-1)]
}  
max.time <- max(Nendtime)
sim.time <- seq(0,max.time,1)-(max.time-60)
duration <- Nendtime
                  
#write max precip values for each plot
max.precip <- as.data.frame(cbind(case.names,Nmaxprecip))
write.csv(max.precip,file=paste(execdir.root,"max.precip.csv",sep=""))

#now we need precip synced to sim.time assuming constant rainfall (Nmaxprecip[])
# estimate precip every minute, ending at minute 60 and estimating for negative time(before flow is observed)
rainfall.mat <- matrix(data = NA, nrow = Nprecip, ncol = length(sim.time))
for(i in 1:Nprecip){
  start.rain <- 60-Nendtime[i]
  end.rain <- Nendtime[i]-(Nendtime[i]-60)
  bump.rain <- (Nendtime[i]-(Ntimesteps[i]-3)*10) #offset from 0
  rain.time.seq <- seq(start.rain,end.rain,1)
  rain.cumulative.seq <- seq(0,Nmaxprecip[i],length.out=length(rain.time.seq))
  #write rain time series to csv
  rain.df <- as.data.frame(cbind(rain.time.seq,rain.cumulative.seq))
  write.csv(rain.df,file=paste(execdir.root,"rain.",case.names[i],".csv",sep=""))
  #write to rainfall.mat
  start.write <- length(sim.time)-length(rain.time.seq)+1
  end.write <- length(sim.time)
  rainfall.mat[i,start.write:end.write] <- rain.cumulative.seq
}
#write entire rain time series matrix to csv
rainfall.df <- as.data.frame(cbind(case.names,rainfall.mat))
dim(rainfall.df)
colnames(rainfall.df) <- c("plot",paste("t=",-218:60,sep=""))
write.csv(rainfall.df,file=paste(execdir.root,"rainfall.all.ts.csv",sep=""))

#extract observed flow data from plots and sync to sim.time
#read observed flow
obs.flow <- read.table(paste(keewokdir,"flow_obs.in",sep=""),sep="\t")
dim(obs.flow)
# units are mm^3/mm^2- cumulative
# these are observations for 5:60 by 5
# find what columns these map to in sim.time
obs.sim.time <- seq(224,279,5)
sim.time[obs.sim.time]
obs.flow.mat <- NULL
obs.flow.mat <- matrix(data=NA,nrow=144,ncol=length(sim.time))
for(i in 1:144){
    obs.flow.mat[i,obs.sim.time] <- as.vector(obs.flow[i,3:14],mode="numeric")
}
#d403 is not in flow file!, something else needs to be deleted- A104, the third row
if(dim(obs.flow.mat)[[1]]==144){obs.flow.mat <- obs.flow.mat[-3,]}
dim(obs.flow.mat)
rownames(obs.flow.mat) <- case.names
colnames(obs.flow.mat) <- c(paste("t=",-218:60,sep=""))
dim(obs.flow.mat) # this is the real observation matrix for cumulative flow in mm3/mm2
#write entire flow time series matrix to csv
obs.flow.df <- as.data.frame(cbind(case.names,obs.flow.mat))
colnames(obs.flow.df) <- c("plot",paste("t=",rain.time.seq,sep=""))  
write.csv(obs.flow.df,file=paste(execdir.root,"obsflow.all.ts.csv",sep=""))
  

#plot cumulative rainfall
#rainfall.pdf <- paste(rundir,"all.rainfall.pdf",sep="")
#par(mfrow=c(2,2))
#pdf(file=rainfall.pdf,width=8.5, height=10.5, bg="white")
#
#dev.off()

print(paste("Finished 02extractRainfallRunoff.R at",date()))




