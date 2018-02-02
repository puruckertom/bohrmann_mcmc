# Hello, world program using doMPI package

# Load the doMPI and Rmpi packages

library("Rmpi")
library("doMPI")

# Initialize the cluster

cl <- startMPIcluster()
registerDoMPI(cl)

# Read the directory names from a file

execdirs <- scan("execdirs","");

results <- foreach(dir=execdirs) %dopar% {
paste("I am",mpi.comm.rank(),"and I see",dir)
}

print(results)

# Tell all slaves to close down, and exit the program

mpi.close.Rslaves(dellog=FALSE)
mpi.quit()
