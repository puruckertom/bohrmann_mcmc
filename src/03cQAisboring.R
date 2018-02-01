#the 143 without D403
length(case.names)
case.names

par.case.names <- paste(par.in$study,par.in$id,sep="")   

length(rownames(obs.flow.mat))
#row.names(obs.flow.arr)


precip.files

colnames(training.sites.parameters)
rownames(training.sites.parameters)

paste(training.sites.parameters$study,training.sites.parameters$id,sep="")

cbind(case.names, par.case.names,rownames(obs.flow.mat),precip.files)
