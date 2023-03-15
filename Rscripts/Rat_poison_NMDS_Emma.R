#load libraries
library(dplyr)
library(vegan)


abacus.dat<-read.csv('Data/ABACUS_plobata_RatPoisionExp.csv', header=T, row.names=1)

#subset NSAF columns
adjnsaf<-select(abacus.dat, contains('ADJNSAF'))

#keep only proteins that have at least 2 unique spectral counts across replicates
nsaf.uniq<-cbind(adjnsaf, abacus.dat$ALL_NUMPEPSUNIQ)
twopeps<-subset(nsaf.uniq, select=X2021_JAN_13_BARKMAN_PLOBATA_RAT_07_ADJNSAF:X2021_JAN_13_BARKMAN_PLOBATA_RAT_41_ADJNSAF, nsaf.uniq[,29]>1)

#keep only coral proteins - remove contaminants
twopeps$protein<-row.names(twopeps)
prot.coral<-subset(twopeps, grepl(paste('m.', collapse="|"), twopeps$protein))
#2901 proteins

#adjust column names
colnames(prot.coral)<-sub('X2021_JAN_13_BARKMAN_PLOBATA_RAT_', "", colnames(prot.coral))
colnames(prot.coral)<-sub('_ADJNSAF', "", colnames(prot.coral))

#add treatment for sample names:
colnames(prot.coral)[c(1,2,28)]<-paste0("FSW1.",1:3)
colnames(prot.coral)[c(6,18)]<-paste0("FSW3.",1:2)
colnames(prot.coral)[c(5,17)]<-paste0("HBP1.",1:2)
colnames(prot.coral)[c(27)]<-paste0("HBP3.",1)
colnames(prot.coral)[c(3,15)]<-paste0("LBP1.",1:2)
colnames(prot.coral)[c(7,19)]<-paste0("LBP3.",1:2)
colnames(prot.coral)[c(4,16)]<-paste0("LIP1.",1:2)
colnames(prot.coral)[c(10,22)]<-paste0("LIP3.",1:2)
colnames(prot.coral)[c(14,26)]<-paste0("MBP1.",1:2)
colnames(prot.coral)[c(8,20)]<-paste0("MBP3.",1:2)
colnames(prot.coral)[c(11,23)]<-paste0("MIP1.",1:2)
colnames(prot.coral)[c(9,21)]<-paste0("MIP3.",1:2)

#remove symb proteins
prot.coral2<-prot.coral[!grepl("symb",prot.coral$protein),]


#perform NMDS
coral.t<-t(prot.coral[,1:28])
coral.tra<-(coral.t+1)
#coral.tra<-data.trans(coral.tra, method='log', plot=F)

nmds.coral<-metaMDS(coral.tra, distance='bray', k=2, trymax=100, autotransform=F)

ordiplot(nmds.coral, choices=c(1,2), type='text', display='sites', cex=0.5)


#NMDS2
coral.t2<-t(prot.coral2[,1:28])
#coral.tra<-(coral.t+1)
nmds.coral2<-metaMDS(coral.t2, distance='bray', k=2, trymax=100, autotransform=F)
ordiplot(nmds.coral2, choices=c(1,2), type='text', display='sites', cex=0.5)



#average technical replicates
FSW1<-cbind(prot.coral[,1], prot.coral[,2], prot.coral[,28])
FSW3<-cbind(prot.coral[,6], prot.coral[,18])
HBP1<-cbind(prot.coral[,5], prot.coral[,17])
HBP3<-cbind(prot.coral[,27]) #REPLICATE FILE IS MISSING 
HIP1<-cbind(prot.coral[,3], prot.coral[,15])
HIP3<-cbind(prot.coral[,12], prot.coral[,24])
LBP1<-cbind(prot.coral[,13], prot.coral[,25])
LBP3<-cbind(prot.coral[,7], prot.coral[,19])
LIP1<-cbind(prot.coral[,4], prot.coral[,16])
LIP3<-cbind(prot.coral[,10], prot.coral[,22])
MBP1<-cbind(prot.coral[,14], prot.coral[,26])
MBP3<-cbind(prot.coral[,8], prot.coral[,20])
MIP1<-cbind(prot.coral[,11], prot.coral[,23])
MIP3<-cbind(prot.coral[,9], prot.coral[,21])

FSW1.avg<-rowMeans(FSW1)
FSW3.avg<-rowMeans(FSW3)
HBP1.avg<-rowMeans(HBP1)
HBP3.avg<-rowMeans(HBP3)
HIP1.avg<-rowMeans(HIP1)
HIP3.avg<-rowMeans(HIP3)
LBP1.avg<-rowMeans(LBP1)
LBP3.avg<-rowMeans(LBP3)
LIP1.avg<-rowMeans(LIP1)
LIP3.avg<-rowMeans(LIP3)
MBP1.avg<-rowMeans(MBP1)
MBP3.avg<-rowMeans(MBP3)
MIP1.avg<-rowMeans(MIP1)
MIP3.avg<-rowMeans(MIP3)

avg.dat<-data.frame(FSW1.avg, FSW3.avg, HBP1.avg, HBP3.avg, HIP1.avg, HIP3.avg, LBP1.avg, LBP3.avg, LIP1.avg, LIP3.avg, MBP1.avg, MBP3.avg, MIP1.avg, MIP3.avg)
rownames(avg.dat)<-rownames(prot.coral)

#write.csv(avg.dat, "Output/Rat/Rat_position_run1_pooled.csv")

avg.t<-t(avg.dat)
avg.tra<-(avg.t+1)
avg.tra<-data.trans(avg.tra, method='log', plot=F)

nmds.avg<-metaMDS(avg.tra, distance='bray', k=2, trymax=100, autotransform=F)

ordiplot(nmds.avg, choices=c(1,2), type='text', display='sites', cex=0.5)


#FSW1 = filtered seawater = grey76
#FSW3 = filtered seawater = grey76
#LIP = low inert pellet = dodgerblue
#MIP = medium inert pellet = dodgerblue2
#HIP = high inert pellet = dodgerblue4
#LBP = low brodifacoum pellet = coral
#MBP = medium brodifacoum pellet = coral2
#HBP = high brodifacoum pellet = coral4

fig.nmds<-ordiplot(nmds.avg, choices=c(1,2), type='none', display='sites')
points(fig.nmds, 'sites', col='grey55', pch=21, bg=c(rep('grey76',2), rep('coral4',2), rep('dodgerblue4',2), rep('coral',2), rep('dodgerblue',2), rep('coral2',2), rep('dodgerblue2',2)), cex=1.5)
legend(x=-0.09, y=.1, legend=c('filtered SW', 'Low Pellet', 'Med. Pellet', 'High Pellet', 'Low RP', 'Med. RP', 'High RP'), pch=19, col=c('grey76', 'dodgerblue', 'dodgerblue2', 'dodgerblue4', 'coral', 'coral2', 'coral4'))

#ANOSIM
treatment<-c('FSW', 'FSW', 'RP', 'RP', 'IP', 'IP', 'RP', 'RP', 'IP', 'IP', 'RP', 'RP', 'IP', 'IP')
hi.lo.treatment<-c('FSW', 'FSW', 'HRP', 'HRP', 'HIP', 'HIP', 'LRP', 'LRP', 'LIP', 'LIP', 'MRP', 'MRP', 'MIP', 'MIP')
avg.row<-data.stand(avg.t, method='total', margin='row', plot=F)
avg.d<-vegdist(avg.row, 'bray')
treatment.anosim<-anosim(avg.d, grouping=treatment)
summary(treatment.anosim)

ANOSIM statistic R: -0.08387 
      Significance: 0.698 

Permutation: free
Number of permutations: 999

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.162 0.210 0.262 0.328 

Dissimilarity ranks between and within classes:
        0%   25%  50%   75% 100%  N
Between  2 20.75 44.5 63.25   91 60
FSW      1  1.00  1.0  1.00    1  1
IP       3 33.50 62.0 71.00   88 15
RP       5 26.50 48.0 68.50   89 15

hi.lo.anosim<-anosim(avg.d, grouping=hi.lo.treatment)
summary(hi.lo.anosim)

ANOSIM statistic R: 0.466 
      Significance: 0.001 

Permutation: free
Number of permutations: 999

Upper quantiles of permutations (null model):
  90%   95% 97.5%   99% 
0.170 0.228 0.262 0.306 

Dissimilarity ranks between and within classes:
        0%   25%  50%   75% 100%  N
Between  2 26.75 47.5 68.25   91 84
FSW      1  1.00  1.0  1.00    1  1
HIP      7  7.00  7.0  7.00    7  1
HRP     24 24.00 24.0 24.00   24  1
LIP     73 73.00 73.0 73.00   73  1
LRP      6  6.00  6.0  6.00    6  1
MIP      3  3.00  3.0  3.00    3  1
MRP     71 71.00 71.0 71.00   71  1