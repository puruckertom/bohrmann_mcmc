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

previous.parameters <- rep(.2, 14)   #  inits for playing around

min.theta.1 <-0
max.theta.1 <-1
  
min.theta.2 <-0
max.theta.2 <-1

min.theta.3 <-0
max.theta.3 <-1

min.theta.4 <-0
max.theta.4 <-1

min.theta.5 <-0
max.theta.5 <-1

min.theta.6 <-0
max.theta.6 <-1

min.theta.7 <-0
max.theta.7 <-1

min.theta.8 <-0
max.theta.8 <-1

min.theta.9 <-0
max.theta.9 <-1

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
data.structure <- data.structure[-141,]

G.matrix <- matrix(c(0), nrow = 52, ncol = 52)  
           
Z.grand <- matrix(c(0), nrow = (143*12), ncol = 52)
           
for (i in 1:143)  {
  if (as.numeric(data.structure$season[i]) == 1)  {Z.grand[(12*(i-1)+1):(12*(i-1)+12),1]<-1}
  if (as.numeric(data.structure$season[i]) == 2)  {Z.grand[(12*(i-1)+1):(12*(i-1)+12),2]<-1}
  if (as.numeric(data.structure$season[i]) == 3)  {Z.grand[(12*(i-1)+1):(12*(i-1)+12),3]<-1}
  if (as.numeric(data.structure$season[i]) == 4)  {Z.grand[(12*(i-1)+1):(12*(i-1)+12),4]<-1}
           
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]))==c(1,1))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),5]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]))==c(1,2))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),6]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]))==c(1,3))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),7]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]))==c(2,1))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),8]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]))==c(2,2))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),9]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]))==c(2,3))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),10]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]))==c(3,1))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),11]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]))==c(3,2))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),12]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]))==c(3,3))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),13]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]))==c(4,1))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),14]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]))==c(4,2))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),15]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]))==c(4,3))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),16]<-1}
           
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(1,1,1))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),17]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(1,1,2))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),18]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(1,1,3))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),19]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(1,2,1))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),20]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(1,2,2))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),21]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(1,2,3))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),22]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(1,3,1))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),23]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(1,3,2))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),24]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(1,3,3))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),25]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(2,1,1))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),26]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(2,1,2))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),27]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(2,1,3))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),28]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(2,2,1))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),29]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(2,2,2))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),30]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(2,2,3))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),31]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(2,3,1))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),32]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(2,3,2))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),33]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(2,3,3))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),34]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(3,1,1))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),35]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(3,1,2))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),36]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(3,1,3))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),37]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(3,2,1))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),38]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(3,2,2))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),39]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(3,2,3))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),40]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(3,3,1))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),41]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(3,3,2))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),42]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(3,3,3))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),43]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(4,1,1))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),44]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(4,1,2))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),45]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(4,1,3))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),46]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(4,2,1))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),47]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(4,2,2))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),48]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(4,2,3))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),49]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(4,3,1))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),50]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(4,3,2))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),51]<-1}
  if (all(c(as.numeric(data.structure$season[i]), as.numeric(data.structure$week[i]), as.numeric(data.structure$day[i]))==c(4,3,3))) {Z.grand[(12*(i-1)+1):(12*(i-1)+12),52]<-1}
}         
           
##  This is the grand Z matrix - when only fitting with a subset of the data, we only need those rows of Z.grand
  ##  corresponding to the sites that stay in the fitting process

training.obs <-  NULL

for (j in (which(data.structure[,1]%in%training.sites)))  {
  training.obs <- append(training.obs, seq((12*(j-1)+1),(12*(j-1)+12)))
}

Z.current <- Z.grand[training.obs,]
  
R.matrix <- matrix(c(0), nrow = nrow(Z.current), ncol = nrow(Z.current))  # The R.matrix dims depend on Z.current dims  
R.block2 <- matrix(c(0), nrow = 48, ncol = 48)
  
##  Here we get the observed runoffs in a format to easily compare to the predicted runoffs

obs.flow.new <- obs.flow.mat[,279]

for (i in 1:11)  {
  obs.flow.new <- cbind(obs.flow.mat[,(279-5*i)],obs.flow.new)
}  

#if(is.element(403, training.sites)) {
#  training.indices.total <- append(c(1:140), c(142:144))
#}

#if(is.element(403, training.sites)==FALSE) {
#  training.sites.total <- append(append(append(training.sites, training.sites), training.sites), training.sites)
#}

obs.flow.diff <- matrix(c(0), nrow = 143, ncol = 12)
obs.flow.diff[,1] <- obs.flow.new[,1]
for (i in 2:12)  {
  obs.flow.diff[,i] <- obs.flow.new[,i]-obs.flow.new[,(i-1)]
}

obs.flow.diff.current <- obs.flow.diff[which(data.structure$site %in% training.sites),]

##  Now we draw candidate parameters from the candidate distribution                     
           
varcov.candidate <- matrix(0, nrow = 14, ncol = 14)  #  Here we could introduce the updated covariance matrix for AM
diag(varcov.candidate) <- 1   #  get rid of the too for AM

candidate.parameters <- rtmvnorm(n=1, mean = previous.parameters, sigma = varcov.candidate, lower = mins,
                                 upper = maxes, algorithm = "rejection")
  
## candidate.parameters[1:9] are the thetas for Kineros2, candidate.parameters[10:12] are the 
  ##  three variance components, candidate.parameters[13] is the overall variance parameter and
  ##  candidate.parameters[14] is the AR(1) correlation parameter rho.


  
#############################################################################
###  HERE KINEROS2 WORKS ITS MAGIC AND OUTPUTS "candidate.pred.flow.mat"  ###
#############################################################################



##  Now we need to grab the predictions from the K2 run and do the Metropolis-Hastings comparison
   ##  These are called candidate.pred.flow.mat

#candidate.pred.flow.mat should be of the same dimensions as obs.doff.flow.current but with the predicted DELTA runoffs from K2
                                 
obs.minus.preds.candidate <- stack(obs.flow.diff.current-candidate.pred.flow.mat)
  
diag(G.matrix)[1:4] <- candidate.parameters[10]
diag(G.matrix)[5:12] <- candidate.parameters[11]
diag(G.matrix)[16:52] <- candidate.parameters[12]
  
times <- 1:12
rho <- candidate.parameters[14]
sigma <- candidate.parameters[13]
H <- abs(outer(times, times, "-"))
V <- rho^H
R.block <- (sigma/(1-rho^2))*V  #  These are the blocks for the block-diagonal R.matrix
R.block2[1:12, 1:12] <- R.block
R.block2[13:24, 13:24] <- R.block
R.block2[25:36, 25:36] <- R.block
R.block2[37:48, 37:48] <- R.block
if(is.element(403, training.sites))
  {for (j in 1:(length(training.sites)-1)) {
     R.matrix[(48*(j-1)+1):(48*(j-1)+48),(48*(j-1)+1):(48*(j-1)+48)] <- R.block2
  }
  R.matrix[(48*(length(training.sites)-1)+1):(48*(length(training.sites)-1)+12), (48*(length(training.sites)-1)+1):(48*(length(training.sites)-1)+12)]  <- R.block
  R.matrix[(48*(length(training.sites)-1)+13):(48*(length(training.sites)-1)+24),(48*(length(training.sites)-1)+13):(48*(length(training.sites)-1)+24)]  <- R.block
  R.matrix[(48*(length(training.sites)-1)+25):(48*(length(training.sites)-1)+36),(48*(length(training.sites)-1)+25):(48*(length(training.sites)-1)+36)]  <- R.block
} 
if(is.element(403, training.sites)==FALSE)  {
for (j in 1:length(training.sites))  {
  R.matrix[(48*(j-1)+1):(48*(j-1)+48),(48*(j-1)+1):(48*(j-1)+48)] <- R.block2
}
}
  
varcov.sigrho.candidate <-  Z.current%*%G.matrix%*%t(Z.current)+R.matrix 
  
##  This acceptance probability calculation requires the FULL M-H calculation because the candidate distribution is a truncated
  ##  normal distribution which is NOT symmetric.  

prob.accept.new.thetas <- min(1, (exp(dmnorm(obs.minus.preds.candidate, mean = rep(0, length(obs.minus.preds.candidate)), 
                                  varcov = varcov.sigrho.candidate, log = TRUE)+log(prior(candidate.parameters[10], 
                                  candidate.parameters[11],candidate.parameters[12],candidate.parameters[13]))
                                      + dtmvnorm(candidate.parameters, mean = previous.parameters, sigma = varcov.candidate, log =T)
                                    - dmnorm(obs.minus.preds.previous, mean = rep(0, length(obs.minus.preds.candidate)), 
                                             varcov = varcov.sigrho.previous, log = TRUE) - log(prior(previous.parameters[10], previous.parameters[11],
                                             previous.parameters[12], previous.parameters[13]))- 
                                               dtmvnorm(previous.parameters, mean = as.numeric(candidate.parameters), sigma = varcov.candidate, log =T))))
  
u <- runif(1)

##  Assign the candidate parameters or the previous parameters as the next MCMC draw as per the above
  ##  calculation.  
           
ifelse(u <= prob.accept.new.thetas, parameters.draws[i,] <- candidate.parameters, parameters.draws[i,]
       <- parameters.draws[i-1,])