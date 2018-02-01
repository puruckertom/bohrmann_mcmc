library(coda)

#Nsim <- 5004
#Simnumber = wherever we are
#thin.coda=10
#burnin <- 4

coda.check <- read.csv(paste(outputdir,"parameters.draws.",Simnumber,".csv",sep=""))[burnin:Simnumber,]
dim(coda.check)
colnames(coda.check)
coda.check.mcmc <- mcmc(coda.check)

#View(coda.check)

dim(coda.check)

#calib.manning
#min.theta.1 <-0.05
#max.theta.1 <-0.8

#calib.cv
#min.theta.2 <-0.01
#max.theta.2 <-50

#calib.relief
#min.theta.3 <-0.01
#max.theta.3 <-400

#calib.spacing
#min.theta.4 <-0.001
#max.theta.4 <-2

#calib.interception
#min.theta.5 <-0.01
#max.theta.5 <-300

#calib.ks
#min.theta.6 <-0.2199
#max.theta.6 <-266.3

#calib.g
#min.theta.7 <-0.1
#max.theta.7 <-500

#calib.rock
#min.theta.8 <-0
#max.theta.8 <-0.2

#calib.distribution
#min.theta.9 <-0.14
#max.theta.9 <-1.43

#mins <- c(rep(c(min.theta.1, min.theta.2, min.theta.3, min.theta.4, min.theta.5, min.theta.6, min.theta.7,
 #         min.theta.8, min.theta.9),4), 0, 0, 0, -1)
#maxes <- c(rep(c(max.theta.1, max.theta.2, max.theta.3, max.theta.4, max.theta.5, max.theta.6, max.theta.7,
  #         max.theta.8, max.theta.9),4), Inf, Inf, Inf, 1)
#maxes.finite <- maxes
#maxes.finite[c(37,38,39)] <- c(max(coda.check.mcmc[,90]),max(coda.check.mcmc[,91]),max(coda.check.mcmc[,92]))

if((Simnumber-burnin)>thin.coda){
  seq.coda <- seq(1,(Simnumber-burnin),thin.coda)
  coda.mcmc <- mcmc(coda.check[seq.coda,],thin=thin.coda)
  #View(coda.mcmc)
  dim(coda.mcmc)
  is.mcmc(coda.mcmc)
  newmins <- c(1,mins,mins,mins,mins)
  newmaxes <- c(1,maxes,maxes,maxes,maxes)
  
  #traces of accepted parameters
  pdf(paste(outputdir,"coda.",Simnumber,".all.traces.pdf",sep=""))
    par(mfrow=c(3,2))
    for(i in 1:ncol(coda.mcmc)){
      if(i>1 & i<38){
        traceplot(coda.mcmc[,i],main=paste("trace",colnames(coda.mcmc)[i]),ylim=c(newmins[i],newmaxes[i]))
        densplot(coda.mcmc[,i],main=paste("density",colnames(coda.mcmc)[i]),xlim=c(newmins[i],newmaxes[i]))
      }else{
        traceplot(coda.mcmc[,i],main=paste("trace",colnames(coda.mcmc)[i]))
        densplot(coda.mcmc[,i],main=paste("density",colnames(coda.mcmc)[i]))      
      }
    }
  dev.off()
  
  rejectionRate <- round(rejectionRate(coda.mcmc)[which(colnames(coda.mcmc)=="Acalib.manning")],digits=3)
  
  #trace of numerator and denominator of MH likelihood
  pdf(paste(outputdir,"coda.",Simnumber,".numer.denom.pdf",sep=""))
    par(mfrow=c(1,1))
      trace.numer <- as.vector(coda.mcmc[,which(colnames(coda.mcmc)=="num1")])
      trace.denom <- as.vector(coda.mcmc[,which(colnames(coda.mcmc)=="den1")])
      trace.accept <- as.vector(coda.mcmc[,which(colnames(coda.mcmc)=="accept")])
      trace.x <- as.vector(coda.mcmc[,which(colnames(coda.mcmc)=="X")])
      trace.max <- max(trace.numer,trace.denom)
      trace.min <- min(trace.numer,trace.denom)
      plot(cbind(trace.x,trace.numer),main="thinned trace of likelihoods",ylab="negative log-likelihood",
           col="red",ylim=c(trace.min,trace.max),type="l",sub=paste("rejection rate =",rejectionRate))
      lines(cbind(trace.x,trace.denom),col="blue")
      points(cbind(trace.x,trace.accept*trace.max),pch=3)
  dev.off()
  
  #trace of both proposed and accepted parameters for the 9 calibrated parameters
  colnames(coda.mcmc)
  accepted.parameter.list <- c(2:10,11:19,20:28,29:37)
  
  if(ncol(coda.mcmc)==93){
      
    pdf(paste(outputdir,"coda.",Simnumber,".seasonal.traces.pdf",sep=""),width=10,height=6)
      par(mfcol=c(3,4))
      for(i in 1:9){
        traceplot(coda.check.mcmc[,(53+i)],main=paste("trace",colnames(coda.check.mcmc)[i+1]), ylim = c(min(coda.check.mcmc[,(53+i)]),max(coda.check.mcmc[,(53+i)])), col = "blue")      
        traceplot(coda.check.mcmc[,(1+i)],col="red",add=TRUE)
        densplot(coda.check.mcmc[,(1+i)],main=paste("density-acc",colnames(coda.check.mcmc)[1+i]), col = "red")
        densplot(coda.check.mcmc[,(53+i)],main=paste("density-prop",colnames(coda.check.mcmc)[53+i]), col = "blue")
        traceplot(coda.check.mcmc[,(62+i)],main=paste("trace",colnames(coda.check.mcmc)[i+10]), ylim = c(min(coda.check.mcmc[,(62+i)]),max(coda.check.mcmc[,(62+i)])), col = "blue")
        traceplot(coda.check.mcmc[,(10+i)],col="red",add=TRUE)
        densplot(coda.check.mcmc[,(10+i)],main=paste("density-acc",colnames(coda.check.mcmc)[10+i]), col = "red")
        densplot(coda.check.mcmc[,(62+i)],main=paste("density-prop",colnames(coda.check.mcmc)[62+i]), col = "blue")     
        traceplot(coda.check.mcmc[,(71+i)],main=paste("trace",colnames(coda.check.mcmc)[i+19]), ylim = c(min(coda.check.mcmc[,(71+i)]),max(coda.check.mcmc[,(71+i)])), col = "blue")
        traceplot(coda.check.mcmc[,(19+i)],col="red",add=TRUE)
        densplot(coda.check.mcmc[,(19+i)],main=paste("density-acc",colnames(coda.check.mcmc)[19+i]), col = "red")
        densplot(coda.check.mcmc[,(71+i)],main=paste("density-prop",colnames(coda.check.mcmc)[71+i]), col = "blue")      
        traceplot(coda.check.mcmc[,(80+i)],main=paste("trace",colnames(coda.check.mcmc)[i+28]), ylim = c(min(coda.check.mcmc[,(80+i)]),max(coda.check.mcmc[,(80+i)])), col = "blue")
        traceplot(coda.check.mcmc[,(28+i)],col="red",add=TRUE)      
        densplot(coda.check.mcmc[,(28+i)],main=paste("density-acc",colnames(coda.check.mcmc)[28+i]), col = "red")
        densplot(coda.check.mcmc[,(80+i)],main=paste("density-prop",colnames(coda.check.mcmc)[80+i]), col = "blue") 
    
        #densplot(coda.check.mcmc[,(1+i)],main=paste("density-acc",colnames(coda.check.mcmc)[1+i]), col = "red")
        #densplot(coda.check.mcmc[,(53+i)],main=paste("density-prop",colnames(coda.check.mcmc)[53+i]), col = "blue")
        #densplot(coda.check.mcmc[,(10+i)],main=paste("density-acc",colnames(coda.check.mcmc)[1+i]), col = "red")
        #densplot(coda.check.mcmc[,(62+i)],main=paste("density-prop",colnames(coda.check.mcmc)[53+i]), col = "blue")      
        #densplot(coda.check.mcmc[,(19+i)],main=paste("density-acc",colnames(coda.check.mcmc)[1+i]), col = "red")
        #densplot(coda.check.mcmc[,(71+i)],main=paste("density-prop",colnames(coda.check.mcmc)[53+i]), col = "blue")
        #densplot(coda.check.mcmc[,(28+i)],main=paste("density-acc",colnames(coda.check.mcmc)[1+i]), col = "red")
        #densplot(coda.check.mcmc[,(80+i)],main=paste("density-prop",colnames(coda.check.mcmc)[53+i]), col = "blue")
        #traceplot(coda.mcmc[,53],main=paste("trace",colnames(coda.mcmc)[i]),ylim=c(newmins[i],newmaxes[i]),col="red")
        #traceplot(coda.mcmc[,i],col="blue",add=TRUE)
        #traceplot(coda.mcmc[,53+9],main=paste("trace",colnames(coda.mcmc)[i+9]),ylim=c(newmins[i],newmaxes[i]),col="red")
        #traceplot(coda.mcmc[,i+9],col="blue",add=TRUE)
        #traceplot(coda.mcmc[,53+18],main=paste("trace",colnames(coda.mcmc)[i+18]),ylim=c(newmins[i],newmaxes[i]),col="red")
        #traceplot(coda.mcmc[,i+18],col="blue",add=TRUE)
        #traceplot(coda.mcmc[,53+27],main=paste("trace",colnames(coda.mcmc)[i+27]),ylim=c(newmins[i],newmaxes[i]),col="red")
        #traceplot(coda.mcmc[,i+27],col="blue",add=TRUE)
        #densplot(coda.mcmc[,i],main=paste("density-acc",colnames(coda.mcmc)[i]),xlim=c(newmins[i],newmaxes[i]),col="blue")
        #densplot(coda.mcmc[,i+9],main=paste("density-acc",colnames(coda.mcmc)[i+9]),xlim=c(newmins[i],newmaxes[i]),col="blue")
        #densplot(coda.mcmc[,i+18],main=paste("density-acc",colnames(coda.mcmc)[i+18]),xlim=c(newmins[i],newmaxes[i]),col="blue")
        #densplot(coda.mcmc[,i+27],main=paste("density-acc",colnames(coda.mcmc)[i+27]),xlim=c(newmins[i],newmaxes[i]),col="blue")
        #densplot(coda.mcmc[,53],main=paste("density-prop",colnames(coda.mcmc)[i]),xlim=c(newmins[i],newmaxes[i]),col="red")
        #densplot(coda.mcmc[,53+9],main=paste("density-prop",colnames(coda.mcmc)[i+9]),xlim=c(newmins[i],newmaxes[i]),col="red")
        #densplot(coda.mcmc[,53+18],main=paste("density-prop",colnames(coda.mcmc)[i+18]),xlim=c(newmins[i],newmaxes[i]),col="red")
        #densplot(coda.mcmc[,53+27],main=paste("density-prop",colnames(coda.mcmc)[i+27]),xlim=c(newmins[i],newmaxes[i]),col="red")
      }
    dev.off()
  }
  
  if(ncol(coda.mcmc)== 95) {
    pdf(paste(outputdir,"coda.",Simnumber,".seasonal.traces.pdf",sep=""),width=10,height=6)
    par(mfcol=c(3,4))
    for(i in 1:9){
      traceplot(coda.check.mcmc[,(54+i)],main=paste("trace",colnames(coda.check.mcmc)[i+1]), ylim = c(min(coda.check.mcmc[,(54+i)]),max(coda.check.mcmc[,(54+i)])), col = "blue")      
      traceplot(coda.check.mcmc[,(1+i)],col="red",add=TRUE)
      densplot(coda.check.mcmc[,(1+i)],main=paste("density-acc",colnames(coda.check.mcmc)[1+i]), col = "red")
      densplot(coda.check.mcmc[,(54+i)],main=paste("density-prop",colnames(coda.check.mcmc)[54+i]), col = "blue")
      traceplot(coda.check.mcmc[,(63+i)],main=paste("trace",colnames(coda.check.mcmc)[i+10]), ylim = c(min(coda.check.mcmc[,(63+i)]),max(coda.check.mcmc[,(63+i)])), col = "blue")
      traceplot(coda.check.mcmc[,(10+i)],col="red",add=TRUE)
      densplot(coda.check.mcmc[,(10+i)],main=paste("density-acc",colnames(coda.check.mcmc)[10+i]), col = "red")
      densplot(coda.check.mcmc[,(63+i)],main=paste("density-prop",colnames(coda.check.mcmc)[63+i]), col = "blue")     
      traceplot(coda.check.mcmc[,(72+i)],main=paste("trace",colnames(coda.check.mcmc)[i+19]), ylim = c(min(coda.check.mcmc[,(72+i)]),max(coda.check.mcmc[,(72+i)])), col = "blue")
      traceplot(coda.check.mcmc[,(19+i)],col="red",add=TRUE)
      densplot(coda.check.mcmc[,(19+i)],main=paste("density-acc",colnames(coda.check.mcmc)[19+i]), col = "red")
      densplot(coda.check.mcmc[,(72+i)],main=paste("density-prop",colnames(coda.check.mcmc)[72+i]), col = "blue")      
      traceplot(coda.check.mcmc[,(81+i)],main=paste("trace",colnames(coda.check.mcmc)[i+28]), ylim = c(min(coda.check.mcmc[,(81+i)]),max(coda.check.mcmc[,(81+i)])), col = "blue")
      traceplot(coda.check.mcmc[,(28+i)],col="red",add=TRUE)      
      densplot(coda.check.mcmc[,(28+i)],main=paste("density-acc",colnames(coda.check.mcmc)[28+i]), col = "red")
      densplot(coda.check.mcmc[,(81+i)],main=paste("density-prop",colnames(coda.check.mcmc)[81+i]), col = "blue") 
      
      #densplot(coda.check.mcmc[,(1+i)],main=paste("density-acc",colnames(coda.check.mcmc)[1+i]), col = "red")
      #densplot(coda.check.mcmc[,(53+i)],main=paste("density-prop",colnames(coda.check.mcmc)[53+i]), col = "blue")
      #densplot(coda.check.mcmc[,(10+i)],main=paste("density-acc",colnames(coda.check.mcmc)[1+i]), col = "red")
      #densplot(coda.check.mcmc[,(62+i)],main=paste("density-prop",colnames(coda.check.mcmc)[53+i]), col = "blue")      
      #densplot(coda.check.mcmc[,(19+i)],main=paste("density-acc",colnames(coda.check.mcmc)[1+i]), col = "red")
      #densplot(coda.check.mcmc[,(71+i)],main=paste("density-prop",colnames(coda.check.mcmc)[53+i]), col = "blue")
      #densplot(coda.check.mcmc[,(28+i)],main=paste("density-acc",colnames(coda.check.mcmc)[1+i]), col = "red")
      #densplot(coda.check.mcmc[,(80+i)],main=paste("density-prop",colnames(coda.check.mcmc)[53+i]), col = "blue")
      #traceplot(coda.mcmc[,53],main=paste("trace",colnames(coda.mcmc)[i]),ylim=c(newmins[i],newmaxes[i]),col="red")
      #traceplot(coda.mcmc[,i],col="blue",add=TRUE)
      #traceplot(coda.mcmc[,53+9],main=paste("trace",colnames(coda.mcmc)[i+9]),ylim=c(newmins[i],newmaxes[i]),col="red")
      #traceplot(coda.mcmc[,i+9],col="blue",add=TRUE)
      #traceplot(coda.mcmc[,53+18],main=paste("trace",colnames(coda.mcmc)[i+18]),ylim=c(newmins[i],newmaxes[i]),col="red")
      #traceplot(coda.mcmc[,i+18],col="blue",add=TRUE)
      #traceplot(coda.mcmc[,53+27],main=paste("trace",colnames(coda.mcmc)[i+27]),ylim=c(newmins[i],newmaxes[i]),col="red")
      #traceplot(coda.mcmc[,i+27],col="blue",add=TRUE)
      #densplot(coda.mcmc[,i],main=paste("density-acc",colnames(coda.mcmc)[i]),xlim=c(newmins[i],newmaxes[i]),col="blue")
      #densplot(coda.mcmc[,i+9],main=paste("density-acc",colnames(coda.mcmc)[i+9]),xlim=c(newmins[i],newmaxes[i]),col="blue")
      #densplot(coda.mcmc[,i+18],main=paste("density-acc",colnames(coda.mcmc)[i+18]),xlim=c(newmins[i],newmaxes[i]),col="blue")
      #densplot(coda.mcmc[,i+27],main=paste("density-acc",colnames(coda.mcmc)[i+27]),xlim=c(newmins[i],newmaxes[i]),col="blue")
      #densplot(coda.mcmc[,53],main=paste("density-prop",colnames(coda.mcmc)[i]),xlim=c(newmins[i],newmaxes[i]),col="red")
      #densplot(coda.mcmc[,53+9],main=paste("density-prop",colnames(coda.mcmc)[i+9]),xlim=c(newmins[i],newmaxes[i]),col="red")
      #densplot(coda.mcmc[,53+18],main=paste("density-prop",colnames(coda.mcmc)[i+18]),xlim=c(newmins[i],newmaxes[i]),col="red")
      #densplot(coda.mcmc[,53+27],main=paste("density-prop",colnames(coda.mcmc)[i+27]),xlim=c(newmins[i],newmaxes[i]),col="red")
    }
    dev.off()
    
  }
  
  
  colnames(coda.mcmc)
  
  if(ncol(coda.mcmc)==93){
  
    pdf(paste(outputdir,"coda.",Simnumber,".var.rho.traces.pdf",sep=""),width=7,height=10)
      par(mfrow=c(4,3))
      for(i in 38:41){
        traceplot(coda.mcmc[,i+52],main=paste("trace",colnames(coda.mcmc)[i]),col="red")
        traceplot(coda.mcmc[,i],col="blue",add=TRUE) 
        densplot(coda.mcmc[,i+52],main=paste("density-prop",colnames(coda.mcmc)[i]),col="red")
        densplot(coda.mcmc[,i],main=paste("density-acc",colnames(coda.mcmc)[i]),col="blue")
      }
    dev.off()
  }
  
  if(ncol(coda.mcmc)==95){
    
    pdf(paste(outputdir,"coda.",Simnumber,".var.rho.traces.pdf",sep=""),width=7,height=10)
    par(mfrow=c(4,3))
    for(i in 37:41){
      traceplot(coda.mcmc[,i+53],main=paste("trace",colnames(coda.mcmc)[i]),col="red")
      traceplot(coda.mcmc[,i],col="blue",add=TRUE) 
      densplot(coda.mcmc[,i+53],main=paste("density-prop",colnames(coda.mcmc)[i]),col="red")
      densplot(coda.mcmc[,i],main=paste("density-acc",colnames(coda.mcmc)[i]),col="blue")
    }
    dev.off()
  }
  
  #residual plot
  seasons.resids.1 <- substr(rownames(training.sites.flow), 1,1)  
  seasons.resids.2 <- rep(seasons.resids.1,12)
  stack.obs.flow.training <- stack(as.data.frame(obs.flow.training))[,1] 
  preds.graph <- -1*(obs.minus.preds.previous- stack.obs.flow.training) 
  pdf(paste(outputdir,"residuals.",Simnumber,".pdf", sep=""))
    plot(x = preds.graph[which(seasons.resids.2=="A")], y = obs.minus.preds.previous[which(seasons.resids.2=="A")],col = "red", pch = 1,
	xlim = c(-.5, max(preds.graph)), ylim = c(min(obs.minus.preds.previous), max(obs.minus.preds.previous)), main = paste("Residuals versus Predictions for run",Simnumber),
      ylab = "Model Residual", xlab = "Predicted Runoff")    
    points(x = preds.graph[which(seasons.resids.2=="B")], y = obs.minus.preds.previous[which(seasons.resids.2=="B")],col = "blue", pch = 1,
	xlim = c(-.5, max(preds.graph)), ylim = c(min(obs.minus.preds.previous), max(obs.minus.preds.previous)))
    points(x = preds.graph[which(seasons.resids.2=="C")], y = obs.minus.preds.previous[which(seasons.resids.2=="C")],col = "green", pch = 1,
	xlim = c(-.5, max(preds.graph)), ylim = c(min(obs.minus.preds.previous), max(obs.minus.preds.previous)))
    points(x = preds.graph[which(seasons.resids.2=="D")], y = obs.minus.preds.previous[which(seasons.resids.2=="D")],col = "black", pch = 1,
	xlim = c(-.5, max(preds.graph)), ylim = c(min(obs.minus.preds.previous), max(obs.minus.preds.previous)))
    legend(x = 0, y = (min(obs.minus.preds.previous)+40), legend = c("Season 1","Season 2","Season 3","Season 4" ), pch = 1, col = c("red", "blue", "green", "black"))
    dev.off()
  pdf(paste(outputdir,"residuals.partialscale.",Simnumber,".pdf", sep=""))
    plot(x = preds.graph[which(seasons.resids.2=="A")], y = obs.minus.preds.previous[which(seasons.resids.2=="A")],col = "red", pch = 1,
      xlim = c(-.5, 2), ylim = c(-2,2), main = paste("Residuals versus Predictions for run - Limited Axes",Simnumber),
      ylab = "Model Residual", xlab = "Predicted Runoff")    
    points(x = preds.graph[which(seasons.resids.2=="B")], y = obs.minus.preds.previous[which(seasons.resids.2=="B")],col = "blue", pch = 1,
      xlim = c(-.5, max(preds.graph)), ylim = c(min(obs.minus.preds.previous), max(obs.minus.preds.previous)))
    points(x = preds.graph[which(seasons.resids.2=="C")], y = obs.minus.preds.previous[which(seasons.resids.2=="C")],col = "green", pch = 1,
      xlim = c(-.5, max(preds.graph)), ylim = c(min(obs.minus.preds.previous), max(obs.minus.preds.previous)))
    points(x = preds.graph[which(seasons.resids.2=="D")], y = obs.minus.preds.previous[which(seasons.resids.2=="D")],col = "black", pch = 1,
      xlim = c(-.5, max(preds.graph)), ylim = c(min(obs.minus.preds.previous), max(obs.minus.preds.previous)))
    legend(x = 0, y = -2, legend = c("Season 1","Season 2","Season 3","Season 4" ), pch = 1, col = c("red", "blue", "green", "black"))
    dev.off()

  #Observed versus Predicted scatter plot
  pdf(paste(outputdir,"observed.versus.predicted.",Simnumber,".pdf", sep=""))
    plot(x = stack.obs.flow.training[which(seasons.resids.2=="A")], y = preds.graph[which(seasons.resids.2=="A")], main = paste("Observed versus Predicted Runoff  - MCMC Iteration",Simnumber),
	xlim = c(-.5, max(obs.minus.preds.previous)), ylim = c(min(preds.graph),max(preds.graph)), col = "red", xlab = "Observed Runoff", ylab = "Predicted Runoff", pch = 1)
    points(x = stack.obs.flow.training[which(seasons.resids.2=="B")], y = preds.graph[which(seasons.resids.2=="B")], col = "blue", pch = 1)
    points(x = stack.obs.flow.training[which(seasons.resids.2=="C")], y = preds.graph[which(seasons.resids.2=="C")], col = "green", pch = 1)
    points(x = stack.obs.flow.training[which(seasons.resids.2=="D")], y = preds.graph[which(seasons.resids.2=="D")], col = "black", pch = 1)
    abline(a = 0, b = 1, col = "red")	
    legend(x = 2, y = max(preds.graph), legend = c("Season 1","Season 2","Season 3","Season 4" ), pch = 1, col = c("red", "blue", "green", "black"))
    dev.off()
}
#end conditional
