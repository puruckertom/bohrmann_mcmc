root.kineros <- "c:\\dropbox\\kineros2\\"

#ex1.par        Example parameter file (simple, hypothetical configuration with
#               one plane contributing laterally to a channel, with erosion and
#               sediment transport)               
#ex1.pre        Example rainfall input file
system.time(
  system(paste(root.kineros,"kineros2_ex1.bat",sep=""))
)

#wg11.par       Parameter file for Walnut Gulch subwatershed 11 with 17 elements,
#               as shown on the Kineros2 home web page (runoff only - no erosion
#               and sediment transport modeling)
#               
#4Aug80.pre     Rainfall input file with 10 spatially distributed rain gages
#               (not the same storm shown on the home web page)
system.time(
  system(paste(root.kineros,"kineros2_ex2.bat",sep=""))
)
