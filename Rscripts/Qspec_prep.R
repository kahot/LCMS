##########
library(dplyr)

spec<-read.csv("Data/MAUI_2018.csv", stringsAsFactors = F,header = T)
dt<-read.table("Data/Alex_ABACUS_output.tsv", stringsAsFactors = F,sep = "\t", header = T)

spec.df<-data.frame(PROTID=spec$PROTID)
d2<-dt[-1,c("PROTID","PROTLEN")]
spec.df<-merge(spec.df, d2, by="PROTID")
spec.df<-spec.df[!grepl("gi",spec.df$PROTID),]

samples<-colnames(spec)[3:ncol(spec)]
samples<-substr(samples, start=1, stop=2)
samples<-unique(samples)

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


##### 

# look at the qspec results based on pooled vs none / coral proteins only vs. corals +symbiodinium
olowalu_list<-list.files("Output/", pattern = "fdr.txt")
summary<-data.frame(test=olowalu_list)
protein_up1<-list()
protein_up2<-list()
for (i in 1:length(olowalu_list)){
        dt<-read.table(paste0("Output/", olowalu_list[i]), stringsAsFactors = F, header=T, sep="\t")
        UPinSite1<-dt[dt$LogFoldChange<=-0.5 & dt$Zstatistic<=-2,]
        UPinSite2<-dt[dt$LogFoldChange>= 0.5 & dt$Zstatistic>=2,]
        
        summary[i,"Up_site1"]<-nrow(UPinSite1)
        summary[i,"Up_site2"]<-nrow(UPinSite2)
        
        protein_up1[[i]]<-UPinSite1$Protein
        names(protein_up1)[i]<-substr(olowalu_list[i],start=8, stop=17)
        protein_up2[[i]]<-UPinSite2$Protein
        names(protein_up2)[i]<-substr(olowalu_list[i],start=8, stop=17)
        
}

names(protein_up1)[1:2]<-c("all","corals")
names(protein_up2)[1:2]<-c("all","corals")

combinations<-t(combn(names(protein_up2),2))

Comb_list<-list()
for (i in 1:6){
        Comb_list[[i]]<-paste0(combinations[i,1],".vs.", combinations[i,2])
}


summary2<-data.frame(no.sample1=matrix(nrow=6, ncol=1))
rownames(summary2)<-Comb_list

#intersect(intersect(a,b),c)
#Reduce(intersect, list(a,b,c))

for (i in 1:6){
        s1<-combinations[i,1]
        s2<-combinations[i,2]
        p1<-protein_up2[[s1]]
        p2<-protein_up2[[s2]]
        
       
        common<-c()
        for (k in 1:length(p1)){
                w<-which(p1[k]==p2)
                if (!is_empty(w)) common<-c(common, p1[k])
        }
        #common<-intersect(p1,p2)
        
        
        summary2$no.sample1[i]<-length(p1)
        summary2$no.sample2[i]<-length(p2)
        summary2$no.common[i]<-length(common)
        summary2$ProtID[i]<-paste(common,sep=",", collapse="")
        
}

write.csv(summary, "Output/Olowalu_test_comparison.csv")
write.csv(summary2, "Output/Olowalu_test_comparison2.csv")


