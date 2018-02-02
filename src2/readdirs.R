library(foreach)
library(doMC)
execdirs <- scan("execdirs","");

registerDoMC()
foreach(dir=execdirs) %dopar% {
   cat(paste(dir,"\n"))
}
