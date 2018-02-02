Sys.info()
#
#  Other print options instead of cat:
#  print -- prints string surrounded by quotes and with leading number
#  print.noquote -- prints string surrounded by quotes
#
cat(paste("Started R test at",date(),"\n"))
if(.Platform$OS.type=="unix"){
  rundir <- path.expand("~/Dropbox/kineros2/")
  cat(paste("Run directory set to ",rundir,"\n"))
  rundir <- path.expand("../kineros2/")
  cat(paste("Run directory set to ",rundir,"\n"))
}
memsize <- as.double(system("freem.csh",intern=TRUE))/1024
cat(paste("Memory size = ",memsize,"GB\n"))
