library("Rmpi")

# Read a list of directory names

execdirs <- scan("execdirs","");

numfile <- function(string) {
   myid = mpi.comm.rank()
   return(paste(myid,string,sep=": "))
}

# Spawn as many slaves as possible

mpi.spawn.Rslaves(needlog=FALSE)

# Farm out the directory names to the slaves

results <- mpi.applyLB(execdirs,numfile)

print(results)

# Tell all slaves to close down, and exit the program

mpi.close.Rslaves(dellog=FALSE)
mpi.quit()
