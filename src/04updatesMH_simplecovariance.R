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


# builds a block matrix whose diagonals are the square matrices provided.
# m1=matrix(runif(10*10),nrow=10,ncol=10)
# m2=matrix(runif(5*5),nrow=5,ncol=5)
# blockMatrix<-blockMatrixDiagonal(m1,m2,m2,m1)
# or
# blockMatrix<-blockMatrixDiagonal(list(m1,m2,m2,m1))
# C.Ladroue

blocks.list <- as.list(rep("holder", length(training.sites)))

blockMatrixDiagonal<-function(matrixList){  
  
  if(is.list(matrixList[[1]])) matrixList<-matrixList[[1]]
  
  dimensions<-sapply(matrixList,FUN=function(x) dim(x)[1])
  finalDimension<-sum(dimensions)
  finalMatrix<-matrix(0,nrow=finalDimension,ncol=finalDimension)
  index<-1
  for(k in 1:length(dimensions)){
    finalMatrix[index:(index+dimensions[k]-1),index:(index+dimensions[k]-1)]<-matrixList[[k]]
    index<-index+dimensions[k]
  }
  finalMatrix
}

#these need to be updated
colnames(training.sites.parameters)
dim(training.sites.parameters)

cols.to.update <- c(4,5,7,8,9,11,12,14,15)
#Keewok's means
#init.9 <- colMeans(training.sites.parameters[,cols.to.update])
# was     c(0.4081869,5.6096374,5.4243243,0.1717542,63.9956673,41.9687832,75.9935234,0.0150757,0.9888738)
# now from Tommy after a 10k burnin
#init.9 <- c(0.5923404,2.227448, 1.305273, 0.4955102,94.28374,  47.81898,  79.95650,  0.1690111,0.1554680)
rows.A <- which(par.in[,2]=="A")
rows.B <- which(par.in[,2]=="B")
rows.C <- which(par.in[,2]=="C")
rows.D <- which(par.in[,2]=="D")
#keewoks overfit initial estimates
#init.36 <- c(colMeans(par.in[rows.A,cols.to.update]),colMeans(par.in[rows.B,cols.to.update]),colMeans(par.in[rows.C,cols.to.update]),colMeans(par.in[rows.D,cols.to.update]))
colnames(training.sites.parameters[,cols.to.update])
#mean of the overfitted estimates for the first 9
#previous.parameters <- rep(.2, 40)   #  inits for playing around
#previous.parameters[1:36] <- rep(init.9,4)   #
#previous.parameters[1:36]<- init.36
#previous.parameters[37:40] <- c(80.94, 15.71, 333.386, 0)
previous.parameters <- as.vector(t(read.csv(paste(inputdir,"previous.parameters.csv",sep=""))))


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
parameters.draws <- matrix(c(0), nrow = Nsim, ncol = 94)

mins <- c(rep(c(min.theta.1, min.theta.2, min.theta.3, min.theta.4, min.theta.5, min.theta.6, min.theta.7,
                min.theta.8, min.theta.9),4), 0, 0, 0, 0, -1)
maxes <- c(rep(c(max.theta.1, max.theta.2, max.theta.3, max.theta.4, max.theta.5, max.theta.6, max.theta.7,
                 max.theta.8, max.theta.9),4), Inf, Inf, Inf, Inf, 1)

##  Now we define (informative) lognormal priors for the variance components.  We suppose they have 
## logmean of 0 and sdlog of prior.sd, which for a large value specifies a nonimformative prior

prior.sd <- 10000

##  Need to specify a function called prior() that calculates the prior density value of a value of
##  sigma.a, sigma.b, sigma.c and rho

prior <- function(a,b,c,d)  prod(plnorm(as.vector(c(a,b,c,d)), meanlog = 0, sdlog = prior.sd))


##  Here we get the observed runoffs in a format to easily compare to the predicted runoffs
# obs.flow.cumulative.arr and sim.flow.cumulative.arr
obs.flow.new <- obs.flow.mat[,279]

for (kk in 1:11)  {
  obs.flow.new <- cbind(obs.flow.mat[,(279-5*kk)],obs.flow.new)
}

obs.flow.diff.temp <- matrix(NA, nrow = nrow(obs.flow.new), ncol = 11)
  
for (kk in 1:nrow(obs.flow.new))  {
  obs.flow.diff.temp[kk,] <- diff(obs.flow.new[kk,])
}
  
obs.flow.diff <- cbind(obs.flow.new[,1],obs.flow.diff.temp)


#obs.flow.new contains the actual observations of rainfall 
#in units of mm3/mm2 for all 143 original cases
#obs.flow.training is what we want to compare the simulated cumulative
#flows against
obs.flow.diff.training <- obs.flow.diff[training.seq,]
dim(obs.flow.diff.training)

obs.flow.training <- matrix(NA, nrow = nrow(obs.flow.diff.training), ncol= ncol(obs.flow.diff.training))

for(y in 1:nrow(obs.flow.diff.training)) {
  obs.flow.training[y,] <- cumsum(obs.flow.diff.training[y,])
}

## Box-Cox transformation
#obs.flow.diff.training <- obs.flow.diff.training^lambda



################################################
################################################
################################################
################################################
################################################
### Tommy MH MCMC sampling

Naccepted=0

parameters.draws[1,1:41] <- as.matrix(previous.parameters)
colnames(parameters.draws) <- c(paste("A",colnames(training.sites.parameters[,cols.to.update]),sep=""),paste("B",colnames(training.sites.parameters[,cols.to.update]),sep=""),
                                paste("C",colnames(training.sites.parameters[,cols.to.update]),sep=""),paste("D",colnames(training.sites.parameters[,cols.to.update]),sep=""),
                                c("varA","varB","varC", "varD", "rho"),c("u","pr_accept","numer","denom","accept"),c("num1","num2","num3","den1","den2","den3","numden"),
                                paste("Aprop",colnames(training.sites.parameters[,cols.to.update]),sep=""),paste("Bprop",colnames(training.sites.parameters[,cols.to.update]),sep=""),
                                paste("Cprop",colnames(training.sites.parameters[,cols.to.update]),sep=""),paste("Dprop",colnames(training.sites.parameters[,cols.to.update]),sep=""),
                                c("prop.varA","prop.varB","prop.varC", "prop.varD", "prop.rho")
                                )

varcov.candidate <- as.matrix(read.csv(paste(inputdir,"varcov.candidate.csv",sep=""), row.names = 1))


colnames(varcov.candidate) <- colnames(parameters.draws)[1:ncol(varcov.candidate)]
rownames(varcov.candidate) <- colnames(parameters.draws)[1:ncol(varcov.candidate)]

var.scale <- .14
#var.scale.init <- .01
varcov.candidate <- varcov.candidate*var.scale

epsilon <- .000001
ident.mat <- matrix(0, nrow = 41, ncol = 41)
diag(ident.mat) <- 1
epsilon.addition <- var.scale*epsilon*ident.mat

#setup matices to receive nseff values
nseff.matrix <- matrix(data=NA,nrow=Nsim+1,ncol=Ntrainingcases)
nseff.matrix.seasonal <- matrix(data=NA,nrow=Nsim+1,ncol=Ntrainingcases)
nseff.matrix.all <- matrix(data=NA,nrow=Nsim+1,ncol=Ntrainingcases)

#6503
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
fail.count = 0

extract.seq <- seq(224,279,5)

###  FIVE MINUTE SUM FUNCTION

five.minute.sum <- function(row) {
  temp.five.sum <- rep(NA, length = 55)
  for (i in 1:55) {
    temp.five.sum[i] <- sum(row[(5*i):(5*(i+1)-1)], na.rm = TRUE) 
  }
  temp.five.sum
}

##  LAMBDA TRANSFORM THE OBSERVATIONS AND MAKE THEM CUMULATIVE

obs.flow.training.lambda <- obs.flow.diff.training^lambda

obs.cumulative.training.lambda <- t(apply(obs.flow.training.lambda, 1, cumsum))


##  Variance/covariance matrix stuff

R.matrix <- matrix(c(0), nrow = nrow(obs.flow.training.lambda)*12, ncol = nrow(obs.flow.training.lambda)*12)  # The R.matrix dims depend on Z.current dims  

times <- 1:12
H <- abs(outer(times, times, "-"))

A.block <- matrix(c(1, rep(0, 11),1,1,rep(0,10), 1,1,1, rep(0,9),rep(1,4),rep(0,8), rep(1,5), rep(0,7),
                    rep(1,6), rep(0,6), rep(1,7), rep(0,5), rep(1,8), 0,0,0,0,rep(1,9),0,0,0,
                    rep(1,10),0,0, rep(1,11), 0, rep(1,12)), nrow = 12, ncol = 12, byrow = T)

A.list <- list()
for(i in 1:length(training.sites)) {
  A.list[[i]] <- A.block
  
}

A <- blockMatrixDiagonal(A.list)
stuck.iter <- 0
stuck.list <- list()

for(iii in 2:Nsim){
  
  print(paste("Starting loop:",iii,date()))
  
  ##  Now we draw candidate parameters from the candidate distribution                     
  
  #varcov.candidate <- matrix(0, nrow = 40, ncol = 40)  #  Here we could introduce the updated covariance matrix for AM
  
  #var.scale <- .002    # variance for the thetas and rho of jumping distribution - range*var.scale is variance
  
  #diag(varcov.candidate) <- c(rep(c((maxes[1]-mins[1])*var.scale,(maxes[2]-mins[2])*var.scale,(maxes[3]-mins[3])*var.scale,
  #  				(maxes[4]-mins[4])*var.scale,(maxes[5]-mins[5])*var.scale,(maxes[6]-mins[6])*var.scale,
  #					(maxes[7]-mins[7])*var.scale,(maxes[8]-mins[8])*var.scale,(maxes[9]-mins[9])*var.scale),
  #					4),1,1,1,var.scale)   #  get rid of this too for AM
  
  if(iii %in% c(start.amc:stop.amc)){varcov.candidate <- cov(parameters.draws[(seq(from = 3, to = iii-1, by = 1)),1:41])*var.scale + epsilon.addition}
  #CRUDE ADAPTIVE METROPOLIS STEP
  
  candidate.parameters <- rtmvnorm(n=1, mean = as.vector(previous.parameters), sigma = varcov.candidate, lower = mins,
                                   upper = maxes, algorithm = "rejection")
  
  #9 original parameters that we are updating
  #need to be assigned to training.sites.parameters
  update.parameter.matrix <- matrix(data=NA,ncol=9,nrow=Ntrainingcases)
  #we need to know
  #for(rep in 1:Ntrainingcases){update.parameter.matrix[rep,] <- rbind(candidate.parameters[1:9])}
  #loop season
  for(rep in 1:Ntrainingcases){
    if(rep < Ntrainingunique+1){update.parameter.matrix[rep,] <- rbind(candidate.parameters[1:9])
    }else if(rep < 2*Ntrainingunique+1){update.parameter.matrix[rep,] <- rbind(candidate.parameters[10:18])
    }else if(rep < 3*Ntrainingunique+1){update.parameter.matrix[rep,] <- rbind(candidate.parameters[19:27])
    }else{update.parameter.matrix[rep,] <- rbind(candidate.parameters[28:36])      
    }
  }
  
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
  temp.cumulative.arr <- TurnTheCrank(Simnumber)
    
  temp.flow.arr <- array(NA, dim=dim(temp.cumulative.arr))
  temp.cumulative.arr.lambda <- array(NA, dim= dim(temp.cumulative.arr))
  #sim.flow.arr[,] <- temp.flow.arr
  
  for (r in 1:nrow(temp.flow.arr))  {
    temp.flow.arr[r,(279-length(which(is.na(temp.cumulative.arr[r,]) ==FALSE))+1):279] <- c(temp.cumulative.arr[r,which(is.na(temp.flow.arr[r,]) == FALSE)][1],diff(temp.cumulative.arr[r,which(is.na(temp.cumulative.arr[r,]) == FALSE)]))  
    temp.flow.arr[r,(279-length(which(is.na(temp.cumulative.arr[r,]) ==FALSE))+1)] <- temp.cumulative.arr[r,((which(is.na(temp.flow.arr[r,]) == FALSE))[1]-1)]
  }  
    
  ###  OK here I need to put flows on the same time step as the data (every five minutes) 
      ###  before I do the transformation, or else it just ain't gonna make sense.  Hence, 
      ###  there should be five minute summed predicted flows which correspond to the data
      ###  which THEN can be transformed.  Otherwise I'm transforming the one-minute flows
      ###  and then later summing those, which is wrong.  DELICATE.  IS THERE AND APPLY FUNCTION
      ###  THAT WILL DO THIS? 
  
  
  five.minute.flows <- t(apply(temp.flow.arr,1,five.minute.sum)) 
  
  
    
    #temp.cumulative.arr.lambda[r, which(is.na(temp.flow.arr[r,])==FALSE)[1]:279]<- cumsum((temp.flow.arr[r, which(is.na(temp.flow.arr[r,])==FALSE)[1]:279])^lambda)
  
  preds.flow.lambda <- five.minute.flows^lambda
  
  any.fails.count <- length(which(is.na(preds.flow.lambda)))
  if(any.fails.count>0){fail.count=fail.count + any.fails.count}
  
  preds.cumulative.lambda <- t(apply(preds.flow.lambda, 1, cumsum))[,44:55]
  
  ###  So we have our correctly formatted, lambda transformed predictions along with observations  
     ###  Now we have to assemble the candidate variance/covariance matrix based on the four
     ###  candidate variance parameters and the candidate rho parameter
  
  rho <- candidate.parameters[41]
  varA <- candidate.parameters[37]
  varB <- candidate.parameters[38]
  varC <- candidate.parameters[39]
  varD <- candidate.parameters[40]
  V <- rho^H
  
  R.block.A <- A.block%*%((varA/(1-rho^2))*V)%*%t(A.block)
  R.block.B <- A.block%*%((varB/(1-rho^2))*V)%*%t(A.block)
  R.block.C <- A.block%*%((varC/(1-rho^2))*V)%*%t(A.block)
  R.block.D <- A.block%*%((varD/(1-rho^2))*V)%*%t(A.block)
  
 
  if(is.element(403, training.sites.unique)){
    for(i in 1:((length(training.sites)+1)/4)) {
      blocks.list[[i]] <- R.block.A
      }
    for(i in ((length(training.sites)+1)/4+1):((length(training.sites)+1)/2)) {
      blocks.list[[i]] <- R.block.B
      }
    for(i in ((length(training.sites)+1)/2):((length(training.sites)+1)*.75)) {
      blocks.list[[i]] <- R.block.C
    }
    for(i in ((length(training.sites)+1)*.75):length(training.sites))  {
      blocks.list[[i]] <- R.block.D
    }
    
  }else{
    for(i in 1:((length(training.sites))/4)) {
      blocks.list[[i]] <- R.block.A
    }
    for(i in ((length(training.sites))/4+1):((length(training.sites))/2)) {
      blocks.list[[i]] <- R.block.B
    }
    for(i in ((length(training.sites))/2):((length(training.sites))*.75)) {
      blocks.list[[i]] <- R.block.C
    }
    for(i in ((length(training.sites))*.75):length(training.sites))  {
      blocks.list[[i]] <- R.block.D
    }  
      
  }
    
    potential.varcov.matrix.cumulative.temp <- blockMatrixDiagonal(blocks.list)
    
    potential.varcov.matrix.cumulative <- (potential.varcov.matrix.cumulative.temp 
                                           + t(potential.varcov.matrix.cumulative.temp))/2
    
    diag(potential.varcov.matrix.cumulative) <- diag(potential.varcov.matrix.cumulative.temp)    
  
  obs.minus.preds.candidate <- as.vector(t(obs.cumulative.training.lambda)-t(preds.cumulative.lambda))
  residual.matrix <- rep(0, length(obs.minus.preds.candidate))
  
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
    previous.pred.flow.mat <- temp.cumulative.arr
    varcov.sigrho.previous <- potential.varcov.matrix.cumulative
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
    numerator.1 <- dmnorm(obs.minus.preds.candidate, mean = as.vector(residual.matrix), varcov = potential.varcov.matrix.cumulative, log = TRUE)
    #prior <- function(a,b,c,d)  prod(plnorm(as.vector(c(a,b,c,d)), meanlog = 0, sdlog = prior.sd)) -- lognormal product
    numerator.2 <- log(prior(candidate.parameters[37], candidate.parameters[38],candidate.parameters[39], candidate.parameters[40]))
    #truncated multivariate normal density for the jumping distribution
    numerator.3 <- dtmvnorm(candidate.parameters, mean = previous.parameters, sigma = varcov.candidate, log =T, lower = mins, upper = maxes)
    numerator <- numerator.1 + numerator.2 + numerator.3
    denominator.1 <- dmnorm(obs.minus.preds.previous, mean = rep(0, length(obs.minus.preds.candidate)), varcov = varcov.sigrho.previous, log = TRUE)
    denominator.2 <- log(prior(previous.parameters[37], previous.parameters[38],previous.parameters[39], previous.parameters[40]))
    denominator.3 <- dtmvnorm(previous.parameters, mean = as.numeric(candidate.parameters), sigma = varcov.candidate, log =T, lower = mins, upper = maxes)
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
  ifelse(u <= prob.accept.new.thetas, parameters.draws[iii,1:41] <- candidate.parameters, parameters.draws[iii,1:41] <- parameters.draws[iii-1,1:41])  
  
  
  previous.parameters <- parameters.draws[iii,1:41]
  
  ##  This keeps track of any long runs that get "stuck" and reverts back to a 
  ##  previous parameter set if the chain does indeed get stuck.
  
  if(iii>100 && all(previous.parameters == parameters.draws[(iii-50), 1:41])){
    previous.parameters <- parameters.draws[(iii-51), 1:41]
    stuck.iter <- stuck.iter+1
    stuck.list[[stuck.iter]] <- list(iii, previous.parameters)
    write.csv(stuck.list, file = paste(outputdir,Nsim,".stuck.list.csv", sep = ""))
  }
  
  #accepted!
  if(u <= prob.accept.new.thetas){
    obs.minus.preds.previous <- obs.minus.preds.candidate
    varcov.sigrho.previous <- potential.varcov.matrix.cumulative
    previous.pred.flow.mat <- temp.flow.arr
    accepted.flow.arr <- temp.cumulative.arr
  }
  
  parameters.draws[iii,42:46] <- c(u,prob.accept.new.thetas,numerator,denominator,accepted)
  parameters.draws[iii,47:53] <- c(numerator.1,numerator.2,numerator.3,denominator.1,denominator.2,denominator.3,num.den)
  parameters.draws[iii,54:94] <- candidate.parameters
  
  #print out pdf of time series and NS
  if(Simnumber %% 1000 == 0 | Simnumber == 10){
    write.csv(parameters.draws,paste(outputdir,"parameters.draws.",Simnumber,".csv",sep=""))
    write.csv(obs.minus.preds.previous,paste(outputdir,"obs.minus.preds.previous.",Simnumber,".csv",sep=""))
    write.csv(obs.flow.training,paste(outputdir,"lambda.obs.flow.training.",Simnumber,".csv",sep=""))
    write.csv(previous.pred.flow.mat,paste(outputdir,"lambda.previous.pred.flow.mat.",Simnumber,".csv",sep=""))
    
    source(paste(rundir,"05codaInspect.R",sep=""))
    
  }
  
  if(any.fails.count>0){print("FAILFAILFAILFAILFAILFAILFAILFAILFAILFAILFAILFAILFAILFAILFAILFAILFAIL")}
  print(paste("Made it through the likelihood/updates for loop",Simnumber))
  print(paste("Made it through loop",Simnumber))
}
