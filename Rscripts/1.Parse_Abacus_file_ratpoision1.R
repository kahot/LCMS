##########
library(dplyr)

#
ab<-read.csv("Data/Alex_ABACUS_output.tsv", stringsAsFactors = F, header = T)

#keep only proteins that have at least 2 unique spectral counts across replicates
ab<-ab[ab$ALL_NUMPEPSUNIQ >1,] #2917

#remove contaminants
ab<-ab[!grepl("CONTAMINANT",ab$DEFLINE),]

# Select spectral count number columns
spec<-ab[,grepl('PROT|NUMSPECSTOT|DEFLINE', colnames(ab))]

#Rename column names by removing unnecessary strings
colnames(spec)<-gsub("_NUMSPECSTOT", "", colnames(spec))

#remove the tag for sample names
sample.name<-colnames(spec)[5]
remove<-substr(sample.name,1,18)

n<-which(colnames(spec)=="ALL")
colnames(spec)[(n+1):ncol(spec)]<-gsub(remove, '',colnames(spec)[(n+1):ncol(spec)])

#remove coral_ 
colnames(spec)[(n+1):ncol(spec)]<-gsub("CORAL_", '',colnames(spec)[(n+1):ncol(spec)])



#without pooling tech reps
write.csv(spec, "Data/Brodifacoum.Alex.Oct02021_non-pooled.csv", row.names = F)

prot.coral<-spec[,5:32]


## Pool technical replicates and use the average (if applied)

spec_pooled<-spec[,1:3]
spec_pooled$FSW1<-rowMeans(cbind(prot.coral[,1], prot.coral[,2], prot.coral[,28]))
spec_pooled$FSW3<-rowMeans(cbind(prot.coral[,6], prot.coral[,18]))
spec_pooled$HBP1<-rowMeans(cbind(prot.coral[,5], prot.coral[,17]))
spec_pooled$HBP3<-rowMeans(cbind(prot.coral[,27])) #REPLICATE FILE IS MISSING 
spec_pooled$HIP1<-rowMeans(cbind(prot.coral[,3], prot.coral[,15]))
spec_pooled$HIP3<-rowMeans(cbind(prot.coral[,12], prot.coral[,24]))
spec_pooled$LBP1<-rowMeans(cbind(prot.coral[,13], prot.coral[,25]))
spec_pooled$LBP3<-rowMeans(cbind(prot.coral[,7], prot.coral[,19]))
spec_pooled$LIP1<-rowMeans(cbind(prot.coral[,4], prot.coral[,16]))
spec_pooled$LIP3<-rowMeans(cbind(prot.coral[,10], prot.coral[,22]))
spec_pooled$MBP1<-rowMeans(cbind(prot.coral[,14], prot.coral[,26]))
spec_pooled$MBP3<-rowMeans(cbind(prot.coral[,8], prot.coral[,20]))
spec_pooled$MIP1<-rowMeans(cbind(prot.coral[,11], prot.coral[,23]))
spec_pooled$MIP3<-rowMeans(cbind(prot.coral[,9], prot.coral[,21]))

#FSW1 = filtered seawater 
#FSW3 = filtered seawater 
#LIP = low inert pellet 
#MIP = medium inert pellet
#HIP = high inert pellet =
#LBP = low brodifacoum pellet 
#MBP = medium brodifacoum pellet 
#HBP = high brodifacoum pellet 


write.csv(spec_pooled, "Data/Brodifacoum.spec.csv", row.names = F)

