# Hello, world program -- gets number of parallel tasks from the
# job scheduler

# Load the R MPI package if it is not already loaded.

if (!is.loaded("mpi_initialize")) {
   library("Rmpi")
}

# Spawn as many slaves as possible

mpi.spawn.Rslaves(needlog=FALSE)

# Tell all slaves to return a message identifying themselves

results <- mpi.remote.exec(paste("I am",mpi.comm.rank(),"of",mpi.comm.size()))

print(results)

# Tell all slaves to close down, and exit the program

mpi.close.Rslaves(dellog=FALSE)
mpi.quit()
