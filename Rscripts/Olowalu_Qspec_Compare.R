
#extract adjusted spectcoutn total values:  

dt<-read.table("Data/Alex_ABACUS_output.tsv", stringsAsFactors = F,sep = "\t", header = T)
ID<-read.csv("Data/alex_sequence_file.csv", stringsAsFactors = F)
ID$No<-gsub("2019_April_12_coral_","",ID$File.Name)

#Select Adjusted SPECTOT column
spec.dt<-dplyr::select(dt, contains("_NUMSPECSADJ"))

#
spec.names<-gsub("X2019_APRIL_12_BACKMAN_CORAL_","", names(spec.dt))
spec.name2<-gsub("_NUMSPECSADJ", "" , spec.names)

sampleNames<-data.frame(ID=rep(NA,length(spec.name2)))
for (i in 1:length(spec.name2)){
        sampleNames$ID[i]<-paste0(ID$Comment[ID$No==spec.name2[i]],".",i)
}

# Select coral protein only
colnames(spec.dt)<- sampleNames$ID
spec.dt<-cbind(dt[,1:2], spec.dt)
spec.dt<-spec.dt[  !grepl("gi", spec.dt$PROTID),]
spec.dt<-spec.dt[!grepl("symb", spec.dt$PROTID),]

samples<-colnames(spec.dt)[3:ncol(spec.dt)]
samples<-substr(samples, start=1, stop=2)
samples<-unique(samples)

for (i in c(1,3,5)){
        dat<-spec.dt %>% dplyr:: select("PROTID","PROTLEN",grep(samples[i], names(spec.dt)), grep(samples[i+1], names(spec.dt)))
        n<-ncol(dat)-2
        colnames(dat)[3:(2+n/2)]<-0
        colnames(dat)[(3+n/2):(2+n)]<-1
        dat2<-dat
        dat2$sum<-rowSums(dat[3:(2+n)])
        dat<-dat[dat2$sum!=0,]
        write.table(dat, paste0("Output/",samples[i],"vs", samples[i+1],"_coral_adj.txt"), quote=F,sep="\t",row.names = F )


        
        }



### Pool Olowalu technical replicates

dat<-spec.dt %>% dplyr:: select("PROTID","PROTLEN",grep(samples[5], names(spec.dt)), grep(samples[6], names(spec.dt)))
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
write.table(dat2, paste0("Output/",samples[5],"vs", samples[6],"_pooled_corals_adj.txt"), quote=F,sep="\t",row.names = F )

##########
# remove proteins with only 1 occurance

olo.files<-list.files("Output/Olowalu/", pattern="corals.txt$")
#1. coral pooled

for (i in 1:length(olo.files)){
        DF<-read.table(paste0("Output/Olowalu/", olo.files[i]), sep="\t", header=T)
        n<-ncol(DF)-2
        DF<-DF[rowSums(DF[3:(2+n)])!=1,]
        colnames(DF)[3:(2+n/2)]<-0
        colnames(DF)[(3+n/2):(2+n)]<-1
        write.table(DF,paste0("Output/Olowalu/2.",olo.files[i]), quote=F,sep="\t",row.names = F)
}

