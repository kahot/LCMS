library(rlang)

# look at the qspec results 

results<-list.files("Output/Tungsten/Nonround_remove1/", pattern = ".fdr")

#attach the protein annotation info
pID<-read.csv("Output/Pacuta_annotations_Uniprot.csv", stringsAsFactors = F) 
for (i in 1:length(results)){
        dt<-read.table(paste0("Output/Tungsten/Nonround_remove1/", results[i]), stringsAsFactors = F, header=T, sep="\t")
        test<-gsub(".txt_qspec_paired_fdr", '',results[i])
        
        re1<-merge(dt,pID[,c("ID","Annotation","KOid","UniprotID")],by.x="Protein" ,by.y="ID",all.x = T)
        write.csv(re1, paste0("Output/Tungsten/Nonround_remove1/",test,".Results.csv"), row.names=T)
}


results<-list.files("Output/Tungsten/Results/", pattern = ".csv")
summary<-data.frame(test=results)
SigProteins<-list()
for (i in 1:length(results)){
        dt<-read.csv(paste0("Output/Tungsten/Results/", results[i]), stringsAsFactors = F)
        test<-gsub(".Results.csv", '',results[i])
        
        UPinTreat<-dt[dt$LogFoldChange<=-0.5 & dt$Zstatistic<=-2,]
        UPinControl<-dt[dt$LogFoldChange>= 0.5 & dt$Zstatistic>=2,]
        
        summary[i,"Up_Control"]<-nrow(UPinControl)
        summary[i,"Up_Treat"]<-nrow(UPinTreat)
        
        UPinTreat2<-dt[dt$LogFoldChange<=-0.5 & dt$Zstatistic<=-2&dt$fdr<=0.01,]
        UPinControl2<-dt[dt$LogFoldChange>= 0.5 & dt$Zstatistic>=2&dt$fdr<=0.01,]
        
        summary[i,"Up_Control_fdr0.01"]<-nrow(UPinControl2)
        summary[i,"Up_Treat_fdr0.01"]<-nrow(UPinTreat2)
        
        
        
        sum<-data.frame()
        sum<-rbind(UPinControl2, UPinTreat2)
        sum$Test<-test
        SigProteins[[i]]<-sum
        names(SigProteins)[i]<-test
}

upProteins<-list()
downProteins<-list()
for (i in 1: length(SigProteins)){
        df<-SigProteins[[i]]
        upProteins[[i]]<-df$Annotation[df$LogFoldChange<=0.05]
        downProteins[[i]]<-df$Annotation[df$LogFoldChange>=0.05]
}

intersect(intersect(upProteins[[1]], upProteins[[2]]),upProteins[[3]])

overlaps<-intersect(upProteins[[1]], upProteins[[3]])

DF<-upProt1[upProt1$Annotation %in% overlaps,]


re2<-list.files("Output/Tungsten/", pattern="Results.csv")
for (i in 1: length(re2)){
        dt<-read.csv(paste0("Output/Tungsten/", re2[i]))
        test<-gsub(".Results.csv", '',re2[i])
        
        UPinControl<-dt[dt$LogFoldChange<=-0.5 & dt$Zstatistic<=-2,]
        UPinTreat<-dt[dt$LogFoldChange>= 0.5 & dt$Zstatistic>=2,]
        
        UPinControl<-UPinControl[nchar(UPinControl$UniprotID)>=6,]
        UPinTreat<-UPinTreat[nchar(UPinTreat$UniprotID)>=6,]
        
        #create iPath input file
        df1<-data.frame(ID=paste0("UNIPROT:",UPinControl$UniprotID))
        df1$width<-"w15"
        df1$color<-"#0040ff"
        df2<-data.frame(ID=paste0("UNIPROT:",UPinTreat$UniprotID))
        df2$width<-"w15"
        df2$color<-"ff6600"
        
        df<-rbind(df1,df2)
        write.table(df, paste0("Output/Tungsten/iPath.QspecResults.",test,".txt"), quote=F, sep=" ",row.names = F,col.names = F)
        

        
}

#Compare round number qspec vs. non-round number results

rd<-read.csv("Output/Tungsten/Norounding/Qspec_Pdam.Control.vs.PD30.100m.Reults.csv", row.names = 1)
nonrd<-read.csv("Output/Tungsten/Norounding/Qspec_Pdam.Control.vs.PD30.100m.Reults.csv",row.names = 1)


rd$sum<-rowSums(rd[,3:6])
nonrd$sum<-rowSums(nonrd[,3:4])
rd1<-rd[rd$sum!=0,]
nonrd1<-nonrd[nonrd$sum!=0,]

common<-merge(rd1[,c(1:8)], nonrd1[,c(1,3:8)], by="Protein")

#intersect(intersect(a,b),c)
#Reduce(intersect, list(a,b,c))






mixed<-df
coral<-df

mixed$analysis<-"mixed"
coral$analysis<-"coralonly"

combined<-merge(mixed,coral, by="ID", all=T)
common<-merge(mixed,coral, by="ID")

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

#######################################
### Rat poison exp
rat_list<-list.files("Output/Rat/", pattern = "fdr.csv")
summary<-data.frame(test=rat_list)
protein_up1<-list()
protein_up2<-list()
for (i in 1:length(rat_list)){
        dt<-read.csv(paste0("Output/Rat/", rat_list[i]), stringsAsFactors = F, header=T)
        UPinSite1<-dt[dt$LogFoldChange<=-0.5 & dt$Zstatistic<=-2,]
        UPinSite2<-dt[dt$LogFoldChange>= 0.5 & dt$Zstatistic>=2,]
        
        summary[i,"Up_site1"]<-nrow(UPinSite1)
        summary[i,"Up_site2"]<-nrow(UPinSite2)
        
        protein_up1[[i]]<-UPinSite1$Protein
        names(protein_up1)[i]<-substr(rat_list[i],start=1, stop=14)
        protein_up2[[i]]<-UPinSite2$Protein
        names(protein_up2)[i]<-substr(rat_list[i],start=1, stop=14)
        
}

