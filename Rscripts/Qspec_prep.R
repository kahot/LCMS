##########
library(dplyr)

spec<-read.csv("Data/MAUI_2018.csv", stringsAsFactors = F,header = T)
dt<-read.table("Data/Alex_ABACUS_output.tsv", stringsAsFactors = F,sep = "\t", header = T)

spec.df<-data.frame(PROTID=spec$PROTID)
d2<-dt[-1,c("PROTID","PROTLEN")]
spec.df<-merge(spec.df, d2, by="PROTID")
spec.df<-spec.df[!grepl("gi",spec.df$PROTID),]
spec.df<-spec.df[grepl("m.",spec.df$PROTID),]

samples<-colnames(spec)[3:ncol(spec)]
samples<-substr(samples, start=1, stop=2)
samples<-unique(samples)
samples<-samples[c(1:4,6,5)]

for (i in c(1,3,5)){
        dat1<-spec %>% dplyr:: select("PROTID",grep(samples[i], names(spec)), grep(samples[i+1], names(spec)))
        dat<-merge(spec.df,dat1, by="PROTID")
        n<-ncol(dat)-2
        colnames(dat)[3:(2+n/2)]<-0
        colnames(dat)[(3+n/2):(2+n)]<-1
        dat2<-dat
        dat2$sum<-rowSums(dat[3:(2+n)])
        dat<-dat[dat2$sum!=0,]
        write.table(dat, paste0("Output/",samples[i],"vs", samples[i+1],"_all.txt"), quote=F,sep="\t",row.names = F )
        coral<-dat[!grepl("symb",dat$PROTID),]
        
        write.table(coral, paste0("Output/",samples[i],"vs", samples[i+1],"_corals.txt"), quote=F,sep="\t",row.names = F ) 
}

## _all.txt includes corals and symbiodinium spec count data.
## _corals.txt only contains coral spec count data.


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

