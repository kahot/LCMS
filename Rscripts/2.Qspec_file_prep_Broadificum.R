## Take the output from parse_abacus_file script, run NMDS, and prepare Qspec input files

##########
library(dplyr)
#Read the run 1  file
specs1<-read.csv("Data/Brodifacoum.spec.csv", stringsAsFactors = F)

#Read the run 2  file
specs2<-read.csv("Data/Brodifacoum.spec.run2.pooled.csv", stringsAsFactors = F)


#Rename the columns for pooling run 1 and 2
colnames(specs1)[4:17]<-c("FSW1","FSW2","HBP1","HBP2","HIP1","HIP2","LBP1","LBP2","LIP1","LIP2","MBP1","MBP2","MIP1","MIP2")
colnames(specs2)[3:11]<-c("FSW3","HBP3","HBP4","LIP3","LIP4","MIP3","MIP4","MBP3","MBP4")

spec_common<-merge(specs1, specs2[,c(1,3:11)], by="PROTID")

spec_merged<-merge(specs1, specs2, by="PROTID", all=T)
spec_merged$DEFLINE<-apply(spec_merged[,c("DEFLINE.x","DEFLINE.y")],1, function(x) {if (is.na(x["DEFLINE.x"])) x["DEFLINE.y"]
                                                                   else x["DEFLINE.x"]})
spec_merged2<-spec_merged[,c("PROTID","PROTLEN","DEFLINE","FSW1","FSW2","FSW3","LIP1","LIP2","LIP3","LIP4",
                             "LBP1","LBP2","MIP1","MIP2","MIP3","MIP4","MBP1","MBP2","MBP3","MBP4",
                             "HIP1","HIP2","HBP1","HBP2","HBP3","HBP4")]
spec_merged2[,4:26][is.na(spec_merged2[,4:26])]<-0

#Add the correct protein lengths coral proteins (don't have the symbio info -so remove symb proteins?)
protlens<-read.table("~/programs/Coral_proteins/P_lobata_annotated.faa.fai", sep = "\t" )
colnames(protlens)[1:2]<-c("PROTID","PROTLEN")
protlens<-protlens[,1:2]
#df<-spec_merged2[1:2]
## add the protein length to missing symb proteins
#df$PROTID2<-df$PROTID
#df$PROTID2<-gsub(".*\\|","", df$PROTID2)

spec_merged2<-spec_merged2[,-2]
spec_merged2<-merge(spec_merged2, protlens[,1:2], by="PROTID",all.x = T)
spec_merged2<-spec_merged2[,c(1,26,2:25)]


#lens$PROTLEN<-apply(lens,1, function(x) {if (is.na(x["PROTLEN.x"])) x["PROTLEN.y"]
#                                                        else x["PROTLEN.x"]})
#lens$PROTLEN<-as.integer(lens$PROTLEN)
#spec_merged2$PROTLEN<-lens$PROTLEN

# add the protein length to missing symb proteins
#lens$symb_id<-lens$PROTID
#lens$symb_id<-gsub(".*\\|","", lens$symb_id)


write.csv(spec_merged2,"Output/Rat/Ratpoison_run1_run2_merged.csv")




#Create qspec input files

#1. remove symbio
spec3<-spec_merged2[!grepl("symbB1.EST", spec_merged2$PROTID),]


#plot NMDS
require(vegan)
require(fields)
require(MASS)
require(fossil)
library(ggplot2)

df<-spec3
df<-df[,-c(1:3)]
df<-t(df)
mds1<-metaMDS(df)
ordiplot(mds1, type="text", display="sites")
mdsDF = as.data.frame(scores(mds1))
mdsDF$Sample<-rownames(mdsDF)
mdsDF$Sample1<-substr(mdsDF$Sample, 1,3)
mdsDF$Sample2<-c(1,1,2,1,1,2,2,1,1,1,1,2,2,1,1,2,2,1,1,1,1,2,2)
mdsDF$Sample2<-as.character(mdsDF$Sample2)

ggplot(mdsDF, aes(x=NMDS1,y=NMDS2,color=Sample2)) + 
        geom_point(size=3.5) + theme_bw() +
        theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
        #scale_shape_manual(values=c(19,17), name="Transplant site",labels=c("Nearshore", "Offshore"))+
        xlab("Axis 1") + ylab("Axis 2")+
        scale_color_manual(values=c("darkblue","lightblue"), labels=c("Run1","Run2"))+
        theme(axis.text=element_text(size=10))+
        theme(legend.title = element_blank())
ggsave("Output/Rat/Run1+Run2_Comparison_AllProteins.pdf", width = 5,height = 3.5)        


#use the common proteins only?
spec_common<-spec_common[,c("PROTID","PROTLEN","DEFLINE","FSW1","FSW2","FSW3","LIP1","LIP2","LIP3","LIP4",
                             "LBP1","LBP2","MIP1","MIP2","MIP3","MIP4","MBP1","MBP2","MBP3","MBP4",
                             "HIP1","HIP2","HBP1","HBP2","HBP3","HBP4")]

#1, with symbio
df<-spec_common

#2. without symbio
df<-spec_common[!grepl("symbB1.EST", spec_common$PROTID),]


df<-df[,-c(1:3)]
df<-t(df)
mds1<-metaMDS(df)
ordiplot(mds1, type="text", display="sites")
mdsDF = as.data.frame(scores(mds1))
mdsDF$Sample<-rownames(mdsDF)
mdsDF$Sample1<-substr(mdsDF$Sample, 1,3)
mdsDF$Sample2<-c(1,1,2,1,1,2,2,1,1,1,1,2,2,1,1,2,2,1,1,1,1,2,2)
mdsDF$Sample2<-as.character(mdsDF$Sample2)

ggplot(mdsDF, aes(x=NMDS1,y=NMDS2,color=Sample2)) + 
        geom_point(size=3.5) + theme_bw() +
        theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
        #scale_shape_manual(values=c(19,17), name="Transplant site",labels=c("Nearshore", "Offshore"))+
        xlab("Axis 1") + ylab("Axis 2")+
        scale_color_manual(values=c("darkblue","lightblue"), labels=c("Run1","Run2"))+
        theme(axis.text=element_text(size=10),legend.title = element_blank())
ggsave("Output/Rat/run1+run2/Run1_Run2_Comparison_CommonProteinsOnly.pdf", width = 5,height = 3.5) 



#samples<-colnames(spec_merged2)[4:26]
#first 3 is the control
samplen<-data.frame(treat=c("FSW","LIP","LBP","MIP","MBP","HIP","HBP"), n=c(3,4,2,4,4,2,4))

#Against control for each treatment that has more than 2 samples
#common proteins (all including symbio)
spec3<-spec_common
spec3<-
for (i in 1:(nrow(samplen)-1)){
        
        treat<-samplen$treat[i+1]
        #grab FSW and the treatment samples
        dat<-spec3 %>% dplyr:: select("PROTID","PROTLEN","FSW1","FSW2","FSW3",grep(treat, names(spec3)))
        n<-ncol(dat)-2
        colnames(dat)[3:5]<-0
        colnames(dat)[6:(samplen$n[i+1]+5)]<-1
        
        #remove the proteins that does not occur
        dat2<-dat
        dat2$sum<-rowSums(dat[3:ncol(dat)])
        dat<-dat[dat2$sum!=0,]
        write.table(dat, paste0("Output/Rat/RAT_merged_commonONly_FSW.vs.", treat,".txt"), quote=F,sep="\t",row.names = F )
}



## Use Run 1 only

specs1<-read.csv("Data/Brodifacoum.spec.csv", stringsAsFactors = F)
#remove symb
specs1<-specs1[!grepl("symbB1.EST", specs1$PROTID),]

#Check the protein length
protlens<-read.table("~/programs/Coral_proteins/P_lobata_annotated.faa.fai", sep = "\t" )
colnames(protlens)[1:2]<-c("PROTID","PROTLEN")
protlens<-protlens[,1:2]

check<-merge(specs1, protlens, by="PROTID",all.x = T)
check$Check<-apply(check, 1, function(x) ifelse(x["PROTLEN.x"]==x["PROTLEN.y"], 0,1))
table(check$Check) #all correct


##LOOK at NMDS

df<-specs1[,-c(1:3)]
df<-t(df)
mds1<-metaMDS(df)
ordiplot(mds1, type="text", display="sites")
mdsDF = as.data.frame(scores(mds1))
mdsDF$Sample<-rownames(mdsDF)
mdsDF$Sample1<-substr(mdsDF$Sample, 1,3)

ggplot(mdsDF, aes(x=NMDS1,y=NMDS2,color=Sample1)) + 
        geom_point(size=3.5) + theme_bw() +
        theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
        xlab("Axis 1") + ylab("Axis 2")+
        #scale_color_manual(values=c("darkblue","lightblue"), labels=c("Run1","Run2"))+
        theme(axis.text=element_text(size=10))
ggsave("Output/Rat/qspec_run2_coral_round/NMDSplot.pdf", width = 5,height = 3.5) 





#Compare all High vs. control
colnames(specs1)
samplen<-data.frame(treat=c("FSW","LIP","LBP","MIP","MBP","HIP","HBP"), n=c(2,2,2,2,2,2,2),
                    treatment=c("Control","Low","Low","Mid","Mid","High","High"))

treatments<-c("Low","Mid","High")

source("Rscripts/round2.R")
for (i in 1:length(treatments)){
        treatment<-treatments[i]
        treat<-samplen$treat[samplen$treatment==treatment]
        #grab FSW and the treatment samples
        dat<-specs1 %>% dplyr:: select("PROTID","PROTLEN","FSW1","FSW3",grep(treat[1], names(specs1)), grep(treat[2], names(specs1)))
        n<-ncol(dat)-2
        colnames(dat)[3:4]<-0
        colnames(dat)[5:ncol(dat)]<-1
 
        ##Round the numbers? It should be integers
        dat[,3:ncol(dat)]<-apply(dat[,3:ncol(dat)], 2, function(x) round2(x,0))
        
        #remove the proteins that does not occur
        dat2<-dat
        dat2$sum<-rowSums(dat[3:ncol(dat)])
        dat<-dat[dat2$sum!=0,]
        write.table(dat, paste0("Output/Rat/qspec_run2_coral_round/Rat_FSW.vs.", treatment,".txt"), quote=F,sep="\t",row.names = F )
}

#Compare all treatments:
treatments<-unique(samplen$treat)
treatments<-treatments[-1]
for (i in 1:length(treatments)){
        treatment<-treatments[i]
        treat<-samplen$treat[samplen$treatment==treatment]
        #grab FSW and the treatment samples
        dat<-specs1 %>% dplyr:: select("PROTID","PROTLEN","FSW1","FSW3",grep(treatment, names(specs1)))
        n<-ncol(dat)-2
        colnames(dat)[3:4]<-0
        colnames(dat)[5:ncol(dat)]<-1
        
        ##Round the numbers? It should be integers
        dat[,3:ncol(dat)]<-apply(dat[,3:ncol(dat)], 2, function(x) round2(x,0))
        
        #remove the proteins that does not occur
        dat2<-dat
        dat2$sum<-rowSums(dat[3:ncol(dat)])
        dat<-dat[dat2$sum!=0,]
        write.table(dat, paste0("Output/Rat/qspec_run1_indivTest/Rat_FSW.vs.", treatment,".txt"), quote=F,sep="\t",row.names = F )
}



#Compare pooled High/Mid vs. Control

dat<-specs1 %>% dplyr:: select("PROTID","PROTLEN","FSW1","FSW3",grep("HBP", names(specs1)), grep("MBP", names(specs1)))
colnames(dat)[3:4]<-0
colnames(dat)[5:8]<-1
dat[,3:ncol(dat)]<-apply(dat[,3:ncol(dat)], 2, function(x) round2(x,0))
#remove the proteins that does not occur
dat2<-dat
dat2$sum<-rowSums(dat[3:ncol(dat)])
dat<-dat[dat2$sum!=0,]
#write.table(dat, paste0("Output/Rat/qspec_Run2_CoralOnly/Rat_coralProteins_FSW.vs.High_Med.txt"), quote=F,sep="\t",row.names = F )
write.table(dat, paste0("Output/Rat/qspec_run2_coral_round/Rat_FSW.vs.High_Med.txt"), quote=F,sep="\t",row.names = F )

#Compare all Mid and High vs. control
dat<-specs1 %>% dplyr:: select("PROTID","PROTLEN","FSW1","FSW3",grep("HBP", names(specs1)), grep("MBP", names(specs1)),grep("HIP", names(specs1)), grep("MIP", names(specs1)))
colnames(dat)[3:4]<-0
colnames(dat)[5:12]<-1
dat[,3:ncol(dat)]<-apply(dat[,3:ncol(dat)], 2, function(x) round2(x,0))
#remove the proteins that does not occur
dat2<-dat
dat2$sum<-rowSums(dat[3:ncol(dat)])
dat<-dat[dat2$sum!=0,]
#write.table(dat, paste0("Output/Rat/qspec_Run2_CoralOnly/Rat_coralProteins_FSW.vs.High_Med.txt"), quote=F,sep="\t",row.names = F )
write.table(dat, paste0("Output/Rat/qspec_run2_coral_round/Rat_FSW.vs.allHigh_Med.txt"), quote=F,sep="\t",row.names = F )



#### # Weighted Metric multi-dimensional scaling ####


df2<-decostand(df,method="hellinger")
wcmd1 = (wcmdscale(vegdist(df2), 2, add = TRUE,eig = TRUE))
Wcmd1 = as.data.frame(scores(wcmd1))

samplenames<-rownames(Wcmd1)
samplenames<-substr(samplenames, start=1, stop=3 )
samples<-substr(samplenames, start=1,stop=1)
for (i in 1:length(samples)){
        if (samples[i]=="F") samples[i]<-"Control"
        if (samples[i]=="H") samples[i]<-"High"
        if (samples[i]=="L") samples[i]<-"Low"
        if (samples[i]=="M") samples[i]<-"Medium"
        
}

treats<-substr(samplenames, start=2,stop=2)
for (i in 1:length(treats)){
        if (treats[i]=="B") treats[i]<-"Broadificum"
        if (treats[i]=="I") treats[i]<-"Innert"
        if (treats[i]=="S") treats[i]<-"Control"
}

ggplot(Wcmd1, aes(x=Dim1,y=Dim2,color=samples,shape=treats)) + 
        geom_point(size=4) + theme_bw() +
        theme(legend.text = element_text(size=14), legend.title = element_text(size=16,face = "bold"))+
        labs(color="Concentration",shape="Treatment")+
        xlab("Axis 1") + ylab("Axis 2")+
        theme(axis.text=element_text(size=10))
#        xlim(-0.2,0.18)+ylim(-0.08, 0.08)
ggsave("./Output/Rat/qspec_run2_coral_round/WMD_plot_all.pdf", height = 4, width = 6)











## common proteins between two runs
df<-spec_common
#df<-spec_common[!grepl("symbB1.EST", spec_common$PROTID),]
df<-df[,-c(1:3)]
df<-t(df)
df2<-decostand(df,method="hellinger")
wcmd1 = (wcmdscale(vegdist(df2), 2, add = TRUE,eig = TRUE))
Wcmd1 = as.data.frame(scores(wcmd1))

samplenames<-rownames(Wcmd1)
samplenames<-substr(samplenames, start=1, stop=3 )
samples<-substr(samplenames, start=1,stop=1)
for (i in 1:length(samples)){
        if (samples[i]=="F") samples[i]<-"Control"
        if (samples[i]=="H") samples[i]<-"High"
        if (samples[i]=="L") samples[i]<-"Low"
        if (samples[i]=="M") samples[i]<-"Medium"
        
}

treats<-substr(samplenames, start=2,stop=2)
for (i in 1:length(treats)){
        if (treats[i]=="B") treats[i]<-"Broadificum"
        if (treats[i]=="I") treats[i]<-"Innert"
        if (treats[i]=="S") treats[i]<-"Control"
}

run<-substr(rownames(Wcmd1), start=4,stop=4)
for (i in 1:length(run)){
        if (run[i]=="1"|run[i]=="2") run[i]<-"Run1"
        if (run[i]=="3"|run[i]=="4") run[i]<-"Run2"
}
run<-c(rep("Run 1", times=14,),rep("Run 2", times=9))
Wcmd1$Run<-run
Wcmd1$treat<-treats
ggplot(Wcmd1, aes(x=Dim1,y=Dim2,color=samples,shape=treats)) + 
        geom_point(size=4) + theme_bw() +
        theme(legend.text = element_text(size=12), legend.title = element_text(size=13,face = "bold"))+
        labs(color="Concentration",shape="Treatment")+
        xlab("Axis 1") + ylab("Axis 2")+
        theme(axis.text=element_text(size=10))


ggplot()+
        geom_point(data=Wcmd1[Wcmd1$Run=="Run2",], aes(x=Dim1,y=Dim2,shape=treat), color=1, size=4.5)+
        geom_point(data=Wcmd1, aes(x=Dim1,y=Dim2,color=samples,shape=treats), size=3.8)+ 
        theme_bw() +
        theme(legend.text = element_text(size=12), legend.title = element_text(size=14,face = "bold"))+
        labs(color="Concentration",shape="Treatment")+
        xlab("Axis 1") + ylab("Axis 2")+
        theme(axis.text=element_text(size=10))
ggsave("./Output/Rat/run1+run2/WMD_plot_both_runs_allProteins.pdf", height = 4, width = 6)


##################
### For Run 2 Only ###

specs1<-read.csv("Data/Brodifacoum.spec.run2.pooled.csv", stringsAsFactors = F)
#remove symb
specs1<-specs1[!grepl("symbB1.EST", specs1$PROTID),]

#attach the protein length
protlens<-read.table("~/programs/Coral_proteins/P_lobata_annotated.faa.fai", sep = "\t" )
colnames(protlens)[1:2]<-c("PROTID","PROTLEN")
protlens<-protlens[,1:2]

specs1<-merge(specs1, protlens, by="PROTID",all.x = T)
specs1<-specs1[, c("PROTID", "DEFLINE","PROTLEN","Control", "HBP1","HBP2","LIP1","LIP2","MIP1",
                    "MIP2","MBP1","MBP2")]
colnames(specs1)[4]<-"FSW"
##LOOK at NMDS

df<-specs1[,-c(1:3)]
df<-t(df)
mds1<-metaMDS(df)
ordiplot(mds1, type="text", display="sites")
mdsDF = as.data.frame(scores(mds1))
mdsDF$Sample<-rownames(mdsDF)
mdsDF$Sample1<-substr(mdsDF$Sample, 1,3)

ggplot(mdsDF, aes(x=NMDS1,y=NMDS2,color=Sample1)) + 
        geom_point(size=3.5) + theme_bw() +
        theme(legend.text = element_text(size=12), legend.title = element_blank())+
        xlab("Axis 1") + ylab("Axis 2")+
        #scale_color_manual(values=c("darkblue","lightblue"), labels=c("Run1","Run2"))+
        theme(axis.text=element_text(size=10))
ggsave("Output/Rat/run2/NMDSplot_run2_allSamples.pdf", width = 5,height = 3.5) 


df2<-decostand(df,method="hellinger")
wcmd1 = (wcmdscale(vegdist(df2), 2, add = TRUE,eig = TRUE))
Wcmd1 = as.data.frame(scores(wcmd1))

samplenames<-rownames(Wcmd1)
samplenames<-substr(samplenames, start=1, stop=3 )
samples<-substr(samplenames, start=1,stop=1)
for (i in 1:length(samples)){
        if (samples[i]=="F") samples[i]<-"Control"
        if (samples[i]=="H") samples[i]<-"High"
        if (samples[i]=="L") samples[i]<-"Low"
        if (samples[i]=="M") samples[i]<-"Medium"
        
}

treats<-substr(samplenames, start=2,stop=2)
for (i in 1:length(treats)){
        if (treats[i]=="B") treats[i]<-"Broadificum"
        if (treats[i]=="I") treats[i]<-"Innert"
        if (treats[i]=="S") treats[i]<-"Control"
}

ggplot(Wcmd1, aes(x=Dim1,y=Dim2,color=samples,shape=treats)) + 
        geom_point(size=4) + theme_bw() +
        theme(legend.text = element_text(size=12), legend.title = element_text(size=14,face = "bold"))+
        labs(color="Concentration",shape="Treatment")+
        xlab("Axis 1") + ylab("Axis 2")+
        theme(axis.text=element_text(size=10))
ggsave("./Output/Rat/run2/WMD_plot_allProteins.pdf", height = 4, width = 6)



#Compare all High vs. control
colnames(specs1)
samplen<-data.frame(treat=c("FSW","LIP","LBP","MIP","MBP","HIP","HBP"), n=c(2,2,2,2,2,2,2),
                    treatment=c("Control","Low","Low","Mid","Mid","High","High"))

treatments<-c("Low","Mid","High")

source("Rscripts/round2.R")
for (i in 1:length(treatments)){
        treatment<-treatments[i]
        treat<-samplen$treat[samplen$treatment==treatment]
        #grab FSW and the treatment samples
        dat<-specs1 %>% dplyr:: select("PROTID","PROTLEN","FSW1","FSW3",grep(treat[1], names(specs1)), grep(treat[2], names(specs1)))
        n<-ncol(dat)-2
        colnames(dat)[3:4]<-0
        colnames(dat)[5:ncol(dat)]<-1
        
        ##Round the numbers? It should be integers
        dat[,3:ncol(dat)]<-apply(dat[,3:ncol(dat)], 2, function(x) round2(x,0))
        
        #remove the proteins that does not occur
        dat2<-dat
        dat2$sum<-rowSums(dat[3:ncol(dat)])
        dat<-dat[dat2$sum!=0,]
        write.table(dat, paste0("Output/Rat/qspec_run2_coral_round/Rat_FSW.vs.", treatment,".txt"), quote=F,sep="\t",row.names = F )
}

#Compare all treatments:
treatments<-unique(samplen$treat)
treatments<-treatments[-1]
for (i in 1:length(treatments)){
        treatment<-treatments[i]
        treat<-samplen$treat[samplen$treatment==treatment]
        #grab FSW and the treatment samples
        dat<-specs1 %>% dplyr:: select("PROTID","PROTLEN","FSW1","FSW3",grep(treatment, names(specs1)))
        n<-ncol(dat)-2
        colnames(dat)[3:4]<-0
        colnames(dat)[5:ncol(dat)]<-1
        
        ##Round the numbers? It should be integers
        dat[,3:ncol(dat)]<-apply(dat[,3:ncol(dat)], 2, function(x) round2(x,0))
        
        #remove the proteins that does not occur
        dat2<-dat
        dat2$sum<-rowSums(dat[3:ncol(dat)])
        dat<-dat[dat2$sum!=0,]
        write.table(dat, paste0("Output/Rat/qspec_run1_indivTest/Rat_FSW.vs.", treatment,".txt"), quote=F,sep="\t",row.names = F )
}


