#########################
#Nsim=10- specified in 00runTheShow

#save the theta
#colnames(training.sites.parameters)
#previous.parameters <- training.sites.parameters
#candidate.parameters
#obs.flow.mat # observations
#previous.pred.flow.mat
#candidate.pred.flow.mat 

##  Metropolis-Hastings and potentially Adaptive Metropolis algorithm for 
  ##  obtaining parameter estimates for Kineros2

## Packages:  
  # install.packages("mnormt")
  # library(mnormt)
  # install.packages("tmvtnorm")
  # library(tmvtnorm)
  
##  Assume we have 9 Kineros2 parameters (theta.1 - theta.9) we want to estimate.  
  ##  Further assume they all have uniform priors.  Then we first specify
  ##  the minima and maxima of these parameters.

#these need to be updated
colnames(training.sites.parameters)
dim(training.sites.parameters)

cols.to.update <- c(4,5,7,8,9,11,12,14,15)
#Keewok's means
#init.9 <- colMeans(training.sites.parameters[,cols.to.update])
# was     c(0.4081869,5.6096374,5.4243243,0.1717542,63.9956673,41.9687832,75.9935234,0.0150757,0.9888738)
# now from Tommy after a 10k burnin
init.9 <- c(0.5923404,2.227448, 1.305273, 0.4955102,94.28374,  47.81898,  79.95650,  0.1690111,0.1554680)
colnames(training.sites.parameters[,cols.to.update])
#mean of the overfitted estimates for the first 9
previous.parameters <- rep(.2, 14)   #  inits for playing around
previous.parameters[1:9] <- init.9   #
previous.parameters[10:14] <- c(39.63, 80.94, 15.71, 333.386, -.324)

#need these from Keewok

#calib.manning
min.theta.1 <-0.05
max.theta.1 <-0.8

#calib.cv
min.theta.2 <-0.01
max.theta.2 <-50

#calib.relief
min.theta.3 <-0.01
max.theta.3 <-400

#calib.spacing
min.theta.4 <-0.001
max.theta.4 <-2

#calib.interception
min.theta.5 <-0.01
max.theta.5 <-300

#calib.ks
min.theta.6 <-0.2199
max.theta.6 <-266.3

#calib.g
min.theta.7 <-0.1
max.theta.7 <-500

#calib.rock
min.theta.8 <-0
max.theta.8 <-0.2

#calib.distribution
min.theta.9 <-0.14
max.theta.9 <-1.43


# Matrix to keep track of parameter draws (aka CODA)
parameters.draws <- matrix(c(0), nrow = Nsim, ncol = 26)


mins <- c(min.theta.1, min.theta.2, min.theta.3, min.theta.4, min.theta.5, min.theta.6, min.theta.7,
          min.theta.8, min.theta.9, 0, 0, 0, 0, -1)
maxes <- c(max.theta.1, max.theta.2, max.theta.3, max.theta.4, max.theta.5, max.theta.6, max.theta.7,
           max.theta.8, max.theta.9, Inf, Inf, Inf, Inf, 1)

##  Now we define (informative) lognormal priors for the variance components.  We suppose they have 
  ## logmean of 0 and sdlog of prior.sd, which for a large value specifies a nonimformative prior
           
prior.sd <- 10000
           
##  Need to specify a function called prior() that calculates the prior density value of a value of
  ##  sigma.a, sigma.b, sigma.c and rho
           
prior <- function(a,b,c,d)  prod(plnorm(as.vector(c(a,b,c,d)), meanlog = 0, sdlog = prior.sd))
           
##  Now we define the various parts of the grand mixed model
           
data.structure <- matrix(c(0), nrow = 144, ncol = 5)
data.structure <- as.data.frame(data.structure)          
names(data.structure) <- c("site", "season", "week", "day", "fecal.source")
data.structure$site <- rep(flow.sites, 4) 
data.structure$season <- c(rep(1,36), rep(2,36), rep(3,36), rep(4,36))
data.structure$week <- rep(weeks, 4)
data.structure$day <- rep(days, 4)
data.structure$fecal.source <- rep(fecal.source, 4)           
data.structure <- data.structure[-131,]

G.matrix <- matrix(c(0), nrow = 52, ncol = 52)  
           
Z.grand <- matrix(c(0), nrow = (143*12), ncol = 52)
           
for (ii in 1:143) {
  if (as.numeric(data.structure$season[ii]) == 1)  {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),1]<-1}
  if (as.numeric(data.structure$season[ii]) == 2)  {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),2]<-1}
  if (as.numeric(data.structure$season[ii]) == 3)  {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),3]<-1}
  if (as.numeric(data.structure$season[ii]) == 4)  {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),4]<-1}
           
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]))==c(1,1))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),5]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]))==c(1,2))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),6]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]))==c(1,3))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),7]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]))==c(2,1))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),8]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]))==c(2,2))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),9]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]))==c(2,3))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),10]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]))==c(3,1))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),11]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]))==c(3,2))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),12]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]))==c(3,3))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),13]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]))==c(4,1))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),14]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]))==c(4,2))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),15]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]))==c(4,3))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),16]<-1}
           
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(1,1,1))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),17]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(1,1,2))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),18]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(1,1,3))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),19]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(1,2,1))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),20]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(1,2,2))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),21]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(1,2,3))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),22]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(1,3,1))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),23]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(1,3,2))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),24]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(1,3,3))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),25]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(2,1,1))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),26]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(2,1,2))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),27]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(2,1,3))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),28]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(2,2,1))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),29]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(2,2,2))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),30]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(2,2,3))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),31]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(2,3,1))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),32]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(2,3,2))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),33]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(2,3,3))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),34]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(3,1,1))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),35]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(3,1,2))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),36]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(3,1,3))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),37]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(3,2,1))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),38]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(3,2,2))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),39]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(3,2,3))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),40]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(3,3,1))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),41]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(3,3,2))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),42]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(3,3,3))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),43]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(4,1,1))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),44]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(4,1,2))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),45]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(4,1,3))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),46]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(4,2,1))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),47]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(4,2,2))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),48]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(4,2,3))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),49]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(4,3,1))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),50]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(4,3,2))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),51]<-1}
  if (all(c(as.numeric(data.structure$season[ii]), as.numeric(data.structure$week[ii]), as.numeric(data.structure$day[ii]))==c(4,3,3))) {Z.grand[(12*(ii-1)+1):(12*(ii-1)+12),52]<-1}
}         

         
##  This is the grand Z matrix - when only fitting with a subset of the data, we only need those rows of Z.grand
  ##  corresponding to the sites that stay in the fitting process




training.obs <-  NULL
#training.seq
for (jj in training.seq)  {
  training.obs <- append(training.obs, seq((12*(jj-1)+1),(12*(jj-1)+12)))
}
Z.current <- matrix(data = NA, nrow = length(training.obs), ncol = ncol(Z.grand))
Z.current <- Z.grand[training.obs,]
  
R.matrix <- matrix(c(0), nrow = nrow(Z.current), ncol = nrow(Z.current))  # The R.matrix dims depend on Z.current dims  
R.block2 <- matrix(c(0), nrow = 48, ncol = 48)
  
##  Here we get the observed runoffs in a format to easily compare to the predicted runoffs
# obs.flow.cumulative.arr and sim.flow.cumulative.arr
obs.flow.new <- obs.flow.mat[,279]

for (kk in 1:11)  {
  obs.flow.new <- cbind(obs.flow.mat[,(279-5*kk)],obs.flow.new)
}
#obs.flow.new contains the actual observations of rainfall 
#in units of mm3/mm2 for all 143 original cases
#obs.flow.training is what we want to compare the simulated cumulative
#flows against
obs.flow.training <- obs.flow.new[training.seq,]
dim(obs.flow.training)

#if(is.element(403, training.sites)) {
#  training.indices.total <- append(c(1:140), c(142:144))
#}
#if(is.element(403, training.sites)==FALSE) {
#  training.sites.total <- append(append(append(training.sites, training.sites), training.sites), training.sites)
#}

# the difference stuff is unnecessary
#obs.flow.diff <- matrix(c(0), nrow = 143, ncol = 12)
#obs.flow.diff[,1] <- obs.flow.new[,1]
#for (kki in 2:12)  {
#  obs.flow.diff[,kki] <- obs.flow.new[,kki]-obs.flow.new[,(kki-1)]
#}
#
#dim(obs.flow.diff)
##obs.flow.diff.current <- obs.flow.diff[which(data.structure$site %in% training.sites),]
#obs.flow.diff.current <- obs.flow.diff[training.seq,]
#dim(obs.flow.diff) #107 12
##convert to cumulative
##based on units of mm/hour
#obs.flow.cumulative <- array(data=NA,dim=c(nrow(obs.flow.diff.current),12))
#for(irow in 1:nrow(obs.flow.cumulative)){
#  cumulative.flow <- 0
#  for(jcol in 1:12){
#    cumulative.flow <- cumulative.flow + 
#    obs.flow.cumulative[irow,jcol] <- 
#  }
#}

################################################
################################################
################################################
################################################
################################################
### Tommy MH MCMC sampling

parameters.draws[1,1:14] <- previous.parameters
colnames(parameters.draws) <- c(colnames(training.sites.parameters[,cols.to.update]),c("var1","var2","var3","var4","rho"),
                                c("u","pr_accept","numer","denom","accept"),c("num1","num2","num3","den1","den2","den3","numden"))
#colnames(parameters.draws)[10:14] <- c("var1","var2","var3","var4","rho")
#colnames(parameters.draws)[15:19] <- c("u","pr_accept","numer","denom","accept")
#colnames(parameters.draws)[20:26] <- c("num1","num2","num3","den1","den2","den3","numden")

Naccepted=0
for(iii in 2:Nsim){
  print(paste("Starting loop:",iii,date()))

  ##  Now we draw candidate parameters from the candidate distribution                     
             
  varcov.candidate <- matrix(0, nrow = 14, ncol = 14)  #  Here we could introduce the updated covariance matrix for AM
  
  var.scale <- .01    # variance for the thetas and rho of jumping distribution - range*var.scale is variance
  
  diag(varcov.candidate) <- c((maxes[1]-mins[1])*var.scale,(maxes[2]-mins[2])*var.scale,(maxes[3]-mins[3])*var.scale,(maxes[4]-mins[4])*var.scale,(maxes[5]-mins[5])*var.scale,(maxes[6]-mins[6])*var.scale,(maxes[7]-mins[7])*var.scale,(maxes[8]-mins[8])*var.scale,(maxes[9]-mins[9])*var.scale,1,1,1,1,var.scale)   #  get rid of this too for AM
  
  candidate.parameters <- rtmvnorm(n=1, mean = previous.parameters, sigma = varcov.candidate, lower = mins,
                                   upper = maxes, algorithm = "rejection")

  
  #9 original parameters that we are updating
  #need to be assigned to training.sites.parameters
  update.parameter.matrix <- matrix(data=NA,ncol=9,nrow=Ntrainingcases)
  for(rep in 1:Ntrainingcases){update.parameter.matrix[rep,] <- rbind(candidate.parameters[1:9])}
  training.sites.parameters[1:Ntrainingcases,cols.to.update] <- update.parameter.matrix #candidate.parameters[1:9]

  #dim(training.sites.parameters)
  #View(training.sites.parameters)  
  ## candidate.parameters[1:9] are the thetas for Kineros2, candidate.parameters[10:12] are the 
    ##  three variance components, candidate.parameters[13] is the overall variance parameter and
    ##  candidate.parameters[14] is the AR(1) correlation parameter rho.
  
  # training.sites.parameters[,funny.calib.seq] <- candidate.parameters[1:9] 

  #for(jjj in 1:Ncases){

    #run simulation for each case
    
      
    #############################################################################
    ###  HERE KINEROS2 WORKS ITS MAGIC AND OUTPUTS "candidate.pred.flow.mat"  ###
    #############################################################################
        #dim(sim.flow.arr) #Nsites=107, Ntimesteps=279, Nsims=10 (currently i)
        Simnumber=iii
        #temp.flow.arr can be diff or cumulative depending on turnthecrank function
        temp.flow.arr <- TurnTheCrank(Simnumber)
        #sim.flow.arr[,] <- temp.flow.arr
        dim(temp.flow.arr)
        extract.seq <- seq(224,279,5)
        #convert simulated flow rates to cumulative flows
        candidate.pred.flow.mat <- temp.flow.arr[,extract.seq]
        print(paste("Made it through loop",Simnumber))
  #}
  

  ##  Now we need to grab the predictions from the K2 run and do the Metropolis-Hastings comparison
     ##  These are called candidate.pred.flow.mat
  
  #####IMPORTANT#########
  #candidate.pred.flow.mat should be of the same dimensions as obs.flow.training but with the predicted DELTA runoffs from K2
  #####IMPORTANT#########  
  dim(obs.flow.training)   
  stack.check <- as.data.frame(obs.flow.training-candidate.pred.flow.mat)
  colnames(stack.check) <- colnames(candidate.pred.flow.mat)
  obs.minus.preds.candidate <- stack(stack.check)
  obs.minus.preds.candidate <- as.matrix(as.numeric(obs.minus.preds.candidate[,1]))
  
  diag(G.matrix)[1:4] <- candidate.parameters[10]
  diag(G.matrix)[5:12] <- candidate.parameters[11]
  diag(G.matrix)[16:52] <- candidate.parameters[12]
    
  times <- 1:12
  rho <- candidate.parameters[14]
  sigma <- candidate.parameters[13]
  H <- abs(outer(times, times, "-"))
  V <- rho^H
  R.block.init <- (sigma/(1-rho^2))*V  #  These are the blocks for the block-diagonal R.matrix
  R.block <- matrix(c(0), nrow = nrow(R.block.init), ncol = ncol(R.block.init))
  for (z in 1:nrow(R.block)) {
    for (w in z:nrow(R.block))  {
      R.block[z,w] <- sum(R.block.init[1:z, 1:w])
    }
  }
 
  R.block <- R.block + t(R.block)
  diag(R.block) <- .5*diag(R.block)
  R.block2 <- matrix(c(0), nrow = 48, ncol = 48)
  R.block2[1:12, 1:12] <- R.block
  R.block2[13:24, 13:24] <- R.block
  R.block2[25:36, 25:36] <- R.block
  R.block2[37:48, 37:48] <- R.block  

  if(is.element(403, training.sites.unique))
    {for (jjjj in 1:(length(training.seq.unique)-1)) {
       R.matrix[(48*(jjjj-1)+1):(48*(jjjj-1)+48),(48*(jjjj-1)+1):(48*(jjjj-1)+48)] <- R.block2
    }
    R.matrix[(48*(length(training.seq.unique)-1)+1):(48*(length(training.seq.unique)-1)+12), (48*(length(training.seq.unique)-1)+1):(48*(length(training.seq.unique)-1)+12)]  <- R.block
    R.matrix[(48*(length(training.seq.unique)-1)+13):(48*(length(training.seq.unique)-1)+24),(48*(length(training.seq.unique)-1)+13):(48*(length(training.seq.unique)-1)+24)]  <- R.block
    R.matrix[(48*(length(training.seq.unique)-1)+25):(48*(length(training.seq.unique)-1)+36),(48*(length(training.seq.unique)-1)+25):(48*(length(training.seq.unique)-1)+36)]  <- R.block
  } 
  if(is.element(403, training.sites.unique)==FALSE)  {
  for (jjjjj in 1:length(training.seq.unique))  {
    R.matrix[(48*(jjjjj-1)+1):(48*(jjjjj-1)+48),(48*(jjjjj-1)+1):(48*(jjjjj-1)+48)] <- R.block2
  }
  }
    
  varcov.sigrho.candidate <-  Z.current%*%G.matrix%*%t(Z.current)+R.matrix 
    
  ##  This acceptance probability calculation requires the FULL M-H calculation because the candidate distribution is a truncated
    ##  normal distribution which is NOT symmetric.  
  #
  residual.matrix <- matrix(data=rep(0, nrow(obs.minus.preds.candidate)),nrow=nrow(obs.minus.preds.candidate),ncol=1)
  obs.minus.preds.candidate <- as.vector(obs.minus.preds.candidate)
  if(iii==2){
    prob.accept.new.thetas=1
    numerator=1
    denominator=1
    numerator.1=1
    numerator.2=1
    numerator.3=1
    denominator.1=1
    denominator.2=1
    denominator.3=1
    num.den=1
    obs.minus.preds.previous <- obs.minus.preds.candidate 
    u=0
  }
  if(iii>2){
    #Tommy check the prob.accept.new.thetas
    #prob.accept.new.thetas <- min(1, (exp(dmnorm(obs.minus.preds.candidate, mean = as.vector(residual.matrix), varcov = varcov.sigrho.candidate, log = TRUE)
    #                                      +log(prior(candidate.parameters[10], candidate.parameters[11],candidate.parameters[12],candidate.parameters[13]))
    #                                      + dtmvnorm(candidate.parameters, mean = previous.parameters, sigma = varcov.candidate, log =T)
    #                                      - dmnorm(obs.minus.preds.previous, mean = rep(0, length(obs.minus.preds.candidate)), varcov = varcov.sigrho.previous, log = TRUE) 
    #                                      - log(prior(previous.parameters[10], previous.parameters[11],previous.parameters[12], previous.parameters[13]))
    #                                      - dtmvnorm(previous.parameters, mean = as.numeric(candidate.parameters), sigma = varcov.candidate, log =T))))
    #candidate.parameters <- rtmvnorm(n=1, mean = previous.parameters, sigma = varcov.candidate, lower = mins,upper = maxes, algorithm = "rejection")
    #probability of theta given y of candidate parameters given data with multivariate normal density 
    numerator.1 <- dmnorm(obs.minus.preds.candidate, mean = as.vector(residual.matrix), varcov = varcov.sigrho.candidate, log = TRUE)
    #prior <- function(a,b,c,d)  prod(plnorm(as.vector(c(a,b,c,d)), meanlog = 0, sdlog = prior.sd)) -- lognormal product
    numerator.2 <- log(prior(candidate.parameters[10], candidate.parameters[11],candidate.parameters[12],candidate.parameters[13]))
    #truncated multivariate normal density for the jumping distribution
    numerator.3 <- dtmvnorm(candidate.parameters, mean = previous.parameters, sigma = varcov.candidate, log =T)
    numerator <- numerator.1 + numerator.2 + numerator.3
    denominator.1 <- dmnorm(obs.minus.preds.previous, mean = rep(0, length(obs.minus.preds.candidate)), varcov = varcov.sigrho.previous, log = TRUE)
    denominator.2 <- log(prior(previous.parameters[10], previous.parameters[11],previous.parameters[12], previous.parameters[13]))
    denominator.3 <- dtmvnorm(previous.parameters, mean = as.numeric(candidate.parameters), sigma = varcov.candidate, log =T)
    denominator <- denominator.1 + denominator.2 + denominator.3
    num.den <- exp(numerator - denominator)
    prob.accept.new.thetas <- min(1, num.den)
  }   
  u <- runif(1)
  ##  Assign the candidate parameters or the previous parameters as the next MCMC draw as per the above
    ##  calculation.  

  ifelse(u <= prob.accept.new.thetas,accepted<-1,accepted<-0)
  Naccepted=Naccepted+accepted
  print(paste("simnumber =",iii))
  print(paste("accept freq =",Naccepted,"/",iii-1,": Percent accepted=",round(Naccepted/(iii-1)*100,digits=2)))
  print(paste("u","prob.accept.new.thetas","numerator","denominator","accepted"))
  print(paste(u,prob.accept.new.thetas,numerator,denominator,accepted))
  
  # parameters.draws keeps track of the MCMC draws (i.e. it's the 'coda')
  ifelse(u <= prob.accept.new.thetas, parameters.draws[iii,1:14] <- candidate.parameters, parameters.draws[iii,1:14] <- parameters.draws[iii-1,1:14])  
  previous.parameters <- parameters.draws[iii,1:14]
  
  if(u <= prob.accept.new.thetas){
    obs.minus.preds.previous <- obs.minus.preds.candidate
    varcov.sigrho.previous <- varcov.sigrho.candidate
  }

  if(u<= prob.accept.new.thetas){
    for(bb in cols.to.update){
      training.sites.parameters[,bb]<-candidate.parameters[which(cols.to.update == bb)]                                  
    }
  }   #  if we accept the new params, this assigns them to all sites in training.sites.parameters           
  
  parameters.draws[iii,15:19] <- c(u,prob.accept.new.thetas,numerator,denominator,accepted)
  parameters.draws[iii,20:26] <- c(numerator.1,numerator.2,numerator.3,denominator.1,denominator.2,denominator.3,num.den)
  #next step
  
  #print out pdf of time series and NS
  
  print(paste("Made it through the likelihood/updates for loop",Simnumber))
}
#######
write.csv(parameters.draws,paste(outputdir,"parameters.draws.",Simnumber,".csv",sep=""))
print(paste("Finished 04updatesMH.R at",date()))
