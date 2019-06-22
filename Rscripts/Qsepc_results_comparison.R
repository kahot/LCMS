library(rlang)

# look at the qspec results based on pooled vs non-pooled, and  coral-proteins only vs. corals +symbiodinium
# for Olowalu samples. 

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

combinations<-t(combn(names(protein_up2[1:5]),2))

Comb_list<-list()
for (i in 1:nrow(combinations)){
        Comb_list[[i]]<-paste0(combinations[i,1],".vs.", combinations[i,2])
}


summary2<-data.frame(no.sample1=matrix(nrow=nrow(combinations), ncol=1))
rownames(summary2)<-Comb_list

#intersect(intersect(a,b),c)
#Reduce(intersect, list(a,b,c))

for (i in 1:nrow(combinations)){
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

# Pooling samples identifies more significantly differentially expressed proteins
# Use coral proteins only, and pooled data for Olowalu.




######
olowalu_list<-list.files("Output/Olowalu/", pattern = ".uniprot.csv")
summary<-data.frame(test=olowalu_list)
protein_up1<-list()
protein_up2<-list()
for (i in 1:length(olowalu_list)){
        dt<-read.csv(paste0("Output/Olowalu/", olowalu_list[i]), stringsAsFactors = F, row.names = 1)
        UPinSite1<-dt[dt$LogFoldChange<=-0.5 & dt$Zstatistic<=-2,]
        UPinSite2<-dt[dt$LogFoldChange>= 0.5 & dt$Zstatistic>=2,]
        
        summary[i,"Up_site1"]<-nrow(UPinSite1)
        summary[i,"Up_site2"]<-nrow(UPinSite2)
        
        protein_up1[[i]]<-UPinSite1$Protein
        names(protein_up1)[i]<-substr(olowalu_list[i],start=10, stop=19)
        protein_up2[[i]]<-UPinSite2$Protein
        names(protein_up2)[i]<-substr(olowalu_list[i],start=10, stop=19)
        
}

combinations<-t(combn(names(protein_up2),2))

Comb_list<-list()
for (i in 1:nrow(combinations)){
        Comb_list[[i]]<-paste0(combinations[i,1],".vs.", combinations[i,2])
}


summary2<-data.frame(no.sample1=matrix(nrow=nrow(combinations), ncol=1))
rownames(summary2)<-Comb_list

#intersect(intersect(a,b),c)
#Reduce(intersect, list(a,b,c))

for (i in 1:nrow(combinations)){
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

write.csv(summary, "Output/Olowalu/Olowalu_test_comparison-removedSingleProtein.csv")
write.csv(summary2, "Output/Olowalu/Olowalu_test_comparison2-removedSingleProtein.csv")

#Don't manually adjust the spec counts for the samples with less technical replicates

