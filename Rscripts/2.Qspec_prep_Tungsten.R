##########
library(dplyr)
spec<-read.csv("Output/Tungsten/Pdam_tungsten_pooled_uniprot2.csv")
#dt<-read.table("Data/Alex_ABACUS_output.tsv", stringsAsFactors = F,sep = "\t", header = T)
colnames(spec)[1]<-"PROTID"
spec.df<-spec[,c("PROTID","PROTLEN")]

#Remove the symbionts 
#spec.df<-spec.df[grepl("Pocillopora_acu",spec.df$PROTID),]

samples<-colnames(spec)[4:(ncol(spec)-1)]
#samples<-substr(samples, start=1, stop=3)
#samples<-unique(samples)

#reorder the sample names
samples<-sort(samples)

#only 30 days have replicates
samples1<-samples[grepl("PD30",samples)]
#remove 100b
samples1<-samples1[samples1!="PD30.100b"]
samples1


for (i in c(1,3,5)){
        dat1<-spec %>% dplyr:: select("PROTID",grep(samples1[i], names(spec)), grep(samples1[i+1], names(spec)))
        #control
        dat2<-spec %>% dplyr:: select("PROTID",grep(samples1[7], names(spec)), grep(samples1[8], names(spec)))
        
        dat<-merge(spec.df,merge(dat2,dat1, by="PROTID"), by="PROTID")
        n<-ncol(dat)-2
        colnames(dat)[3:(2+n/2)]<-0
        colnames(dat)[(3+n/2):(2+n)]<-1
        dat2<-dat
        dat2$sum<-rowSums(dat[3:(2+n)])
        #dat<-dat[dat2$sum!=0,]
        #
        ##Round the numbers? It should be integers
        #dat[,3:6]<-apply(dat[,3:6], 2, function(x) ceiling(x))
        #
        ##remove the proteins with single observation of 1
        #dat2<-dat
        #dat2$sum<-rowSums(dat[3:(2+n)])
        dat<-dat[dat2$sum>1,]
        
        #write.table(dat, paste0("Output/Tungsten/Qspec_Pdam.Control.vs.",samples1[i],".coralproteinsonly.txt"), quote=F,sep="\t",row.names = F )
        write.table(dat, paste0("Output/Tungsten/Nonround_remove1/Qspec_Pdam.Control.vs.",samples1[i],".txt"), quote=F,sep="\t",row.names = F )
}



## _all.txt includes corals and symbiodinium spec count data.
## _corals.txt only contains coral spec count data.

