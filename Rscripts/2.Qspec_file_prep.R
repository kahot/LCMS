##########
library(dplyr)

#Continued from Parse_Abacus_file.R
specs<-read.csv("Data/Brodifacoum.spec.csv", stringsAsFactors = F)

specs<-specs[,-3]
samples<-colnames(specs)[3:ncol(specs)]
samples<-substr(samples, start=1, stop=3)
samples<-unique(samples)
#reorder the sample names
samples<-samples[c(1,5,7,3,4,6,2)]
dir.create("Output/Rat")
#compare FSW vs. different treatments
for (i in 2:length(samples)){
        dat<-spec %>% dplyr:: select("PROTID","PROTLEN","FSW1","FSW3",grep(samples[i], names(spec)))
        n<-ncol(dat)-2
        colnames(dat)[3:(2+n/2)]<-0
        colnames(dat)[(3+n/2):(2+n)]<-1
        dat2<-dat
        dat2$sum<-rowSums(dat[3:(2+n)])
        dat<-dat[dat2$sum!=0,]
        write.table(dat, paste0("Output/Rat/RAT_FSW.vs.", samples[i],".txt"), quote=F,sep="\t",row.names = F )
        #coral<-dat[!grepl("symb",dat$PROTID),]
        
        #write.table(coral, paste0("Output/",samples[i],"vs", samples[i+1],"_corals.txt"), quote=F,sep="\t",row.names = F ) 
}

## _all.txt includes corals and symbiodinium spec count data.
## _corals.txt only contains coral spec count data.


#compare innert vs. active
comb<-data.frame(I="LIP")
for (i in 2:4){
        dat<-specs %>% dplyr:: select("PROTID","PROTLEN",grep(samples[i], names(specs)),
                                     grep(samples[i+3], names(specs)))
        n<-ncol(dat)-2
        colnames(dat)[3:(2+n/2)]<-0
        colnames(dat)[(3+n/2):(2+n)]<-1
        dat2<-dat
        dat2$sum<-rowSums(dat[3:(2+n)])
        dat<-dat[dat2$sum!=0,]
        write.table(dat, paste0("Output/Rat/", samples[i],"vs",samples[i+3],".txt"), quote=F,sep="\t",row.names = F )
        #coral<-dat[!grepl("symb",dat$PROTID),]
        
        #write.table(coral, paste0("Output/",samples[i],"vs", samples[i+1],"_corals.txt"), quote=F,sep="\t",row.names = F ) 
}







###
specs<-read.csv("Data/Brodifacoum.spec_non-pooled.csv", stringsAsFactors = F)

spec<-specs[,-c(1:4) ]

colnames(spec)[c(1,2,28)]<-paste0("FSW1.",1:3)
colnames(spec)[c(6,18)]<-paste0("FSW3.",1:2)
colnames(spec)[c(5,17)]<-paste0("HBP1.",1:2)
colnames(spec)[c(27)]<-paste0("HBP3")
colnames(spec)[c(3,15)]<-paste0("HIP1.",1:2)
colnames(spec)[c(12,24)]<-paste0("HIP3.",1:2)
colnames(spec)[c(13,25)]<-paste0("LBP1.",1:2)
colnames(spec)[c(7,19)]<-paste0("LBP3.",1:2)
colnames(spec)[c(4,16)]<-paste0("LIP1.",1:2)
colnames(spec)[c(10,22)]<-paste0("LIP3.",1:2)
colnames(spec)[c(14,26)]<-paste0("MBP1.",1:2)
colnames(spec)[c(8,20)]<-paste0("MBP3.",1:2)
colnames(spec)[c(11,23)]<-paste0("MIP1.",1:2)
colnames(spec)[c(9,21)]<-paste0("MIP1.",1:2)

write.csv(spec,"Data/Brodifacoum.spec_non-pooled_named.csv")


samples

samples<-colnames(spec)
samples<-substr(samples, start=1, stop=3)
samples<-unique(samples)


spec<-spec[,order(colnames(spec))]
spec<-cbind(specs[,1:2],spec)
#dir.create("Output/Rat")
#compare FSW vs. different treatments
for (i in 2:length(samples)){
        dat<-spec %>% dplyr:: select("PROTID","PROTLEN",grep("FSW", names(spec)),,grep(samples[i], names(spec)))
        colnames(dat)[3:7]<-0
        colnames(dat)[8:ncol(dat)]<-1
        dat2<-dat
        dat2$sum<-rowSums(dat[3:ncol(dat)])
        dat<-dat[dat2$sum!=0,]
        write.table(dat, paste0("Output/Rat/RATnonpooled_FSW.vs.", samples[i],".txt"), quote=F,sep="\t",row.names = F )
        #coral<-dat[!grepl("symb",dat$PROTID),]
        
        #write.table(coral, paste0("Output/",samples[i],"vs", samples[i+1],"_corals.txt"), quote=F,sep="\t",row.names = F ) 
}

#compare innert vs. active (non-pooled)
#order samples
samples<-samples[c(1,2,6,5,4,3,7)]

for (i in 2:4){
        dat<-spec %>% dplyr:: select("PROTID","PROTLEN",grep(samples[i], names(spec)),
                                      grep(samples[i+3], names(spec)))

        colnames(dat)[grep(samples[i], names(dat))]<-0
        colnames(dat)[grep(samples[i+3], names(dat))]<-1
        dat2<-dat
        dat2$sum<-rowSums(dat[3:ncol(dat)])
        dat<-dat[dat2$sum!=0,]
        write.table(dat, paste0("Output/Rat/Nonpooled_", samples[i],"vs",samples[i+3],".txt"), quote=F,sep="\t",row.names = F )
        #coral<-dat[!grepl("symb",dat$PROTID),]
        
        #write.table(coral, paste0("Output/",samples[i],"vs", samples[i+1],"_corals.txt"), quote=F,sep="\t",row.names = F ) 
}


###########
### Pool Olowalu technical replicates

samples2<-colnames(spec)[3:ncol(spec)]
samples2<-substr(samples2, start=1, stop=3)
samples2<-unique(samples2)

for (i in 5){
        dat1<-spec %>% dplyr:: select("PROTID",grep(samples[i], names(spec)), grep(samples[i+1], names(spec)))
        dat<-merge(spec.df,dat1, by="PROTID")
        colnames(dat)[3:18]<-substr(colnames(dat)[3:18], start=1, stop=3)
        cnames<-unique(colnames(dat)[3:18])
        dat2<-dat[,1:2]
        
        for (j in 1:length(cnames)){
                dat2[,(2+j)]<-rowSums(dat[,which(colnames(dat)==cnames[j])])
                colnames(dat2)[(2+j)]<-cnames[j]
        }
        n<-ncol(dat2)-2
        colnames(dat2)[3:(2+n/2)]<-0
        colnames(dat2)[(3+n/2):(2+n)]<-1
        dat3<-dat2
        dat3$sum<-rowSums(dat2[3:(2+n)])
        dat2<-dat2[dat3$sum!=0,]
        write.table(dat2, paste0("Output/",samples[i],"vs", samples[i+1],"_pooled_all.txt"), quote=F,sep="\t",row.names = F )
        coral<-dat2[!grepl("symb",dat2$PROTID),]
        write.table(coral, paste0("Output/",samples[i],"vs", samples[i+1],"_pooled_corals.txt"), quote=F,sep="\t",row.names = F ) 
}

########
# eliminate the rows with all zero before running qspec

M2017<-list.files("Output/", pattern="2017.txt$")
for (i in 1:length(olo.files)){
        DF<-read.table(paste0("Output/", M2017[i]), sep="\t", header=T)
        DF<-DF[!grepl("symb",DF$PROTID),]
        DF<-DF[grepl("m.",DF$PROTID),]
        
        n<-ncol(DF)-2
        DF<-DF[rowSums(DF[3:(2+n)])!=0,]
        colnames(DF)[3:(2+n/2)]<-0
        colnames(DF)[(3+n/2):(2+n)]<-1
        write.table(DF,paste0("Output/cleaned.",M2017[i]), quote=F,sep="\t",row.names = F)
}

#h<-read.table("Output/cleaned.WNvsWO2017.txt", stringsAsFactors = F,sep = "\t", header = T)

