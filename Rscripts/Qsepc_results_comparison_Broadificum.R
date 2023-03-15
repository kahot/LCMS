#Take the qspec output file and create a csv file with protein annotation info.

# Here I also looked at how many proteins were significantly more or less abundant 
# between a pair of treatments to understand which input data and qspec runs were most appropriate to analyze your data. 

library(rlang)

# look at the qspec results 
results<-list.files("Output/Rat/qspec_run1+run2_mergedData/", pattern = ".fdr")

summary<-data.frame(test=gsub(".txt_qspec_paired_fdr",'',results))
for (i in 1:length(results)){
        dt<-read.table(paste0("Output/Rat/qspec_run1+run2_mergedData/", results[i]), stringsAsFactors = F, header=T, sep="\t")
        UPinSite1<-dt[dt$LogFoldChange<=-0.5 & dt$Zstatistic<=-2 & dt$fdr<=0.01,]
        UPinSite2<-dt[dt$LogFoldChange>= 0.5 & dt$Zstatistic>=2& dt$fdr<=0.01,]
        
        summary[i,"Up_control"]<-nrow(UPinSite1)
        summary[i,"Up_treatment"]<-nrow(UPinSite2)
}

summary
#                              test Up_control Up_treatment
#1  RAT_merged_commonONly_FSW.vs.HBP        826            0
#2  RAT_merged_commonONly_FSW.vs.HIP       1081            0
#3  RAT_merged_commonONly_FSW.vs.LBP       1522            0
#4  RAT_merged_commonONly_FSW.vs.LIP        861            0
#5  RAT_merged_commonONly_FSW.vs.MBP       1028            0
#6  RAT_merged_commonONly_FSW.vs.MIP       1003            0
#7             RAT_merged_FSW.vs.HBP        918          902
#8             RAT_merged_FSW.vs.HIP        861          914
#9             RAT_merged_FSW.vs.LBP        883          927
#10            RAT_merged_FSW.vs.LIP        935          915
#11            RAT_merged_FSW.vs.MBP        921          910
#12            RAT_merged_FSW.vs.MIP        935          921

#Do NOT pool the runs


#results<-list.files("Output/Rat/qspec_Run1_CoralOnly/", pattern = ".fdr")

results<-list.files("Output/Rat/qspec_run1_coral_round/", pattern = ".fdr")

summary<-data.frame(test=c("High.Med.All","High.Med.broadi","High","Low","Med"))
for (i in 1:length(results)){
        dt<-read.table(paste0("Output/Rat/qspec_run1_coral_round_run2/", results[i]), stringsAsFactors = F, header=T, sep="\t")
        UPinSite1<-dt[dt$LogFoldChange<=-0.5 & dt$Zstatistic<=-2 & dt$fdr<=0.05,]
        UPinSite2<-dt[dt$LogFoldChange>= 0.5 & dt$Zstatistic>=2& dt$fdr<=0.05,]
        
        summary[i,"Up_control"]<-nrow(UPinSite1)
        summary[i,"Up_treatment"]<-nrow(UPinSite2)
}



#attach the protein annotation info
pID<-read.csv("Data/ProteinID_sorted.csv", stringsAsFactors = F) 
Results<-list()
for (i in 1:length(results)){
        dt<-read.table(paste0("Output/Rat/qspec_run1_coral_round_run2/", results[i]), stringsAsFactors = F, header=T, sep="\t")
        test<-gsub(".txt_qspec_paired_fdr", '',results[i])
        
        re1<-merge(dt,pID,by="Protein",all.x = T)
        Results[[i]]<-re1
        names(Results)[i]<-test
       write.csv(re1, paste0("Output/Rat/qspec_run1_coral_round_run2/",test,".Results.csv"), row.names=T)
}


annotation<-function(x){
        def<-pID$DEFLINE[pID$Protein==x]
        return(def)
}

annotation("m.5611")



#check which proteins were more abundnat in treatments
upProteins<-list()
downProteins<-list()
for (i in 1:length(results)){
        dt<-Results[[i]]
        test<-gsub(".Results.csv", '',names(Results)[i])
        
        UPinControl<-dt[dt$LogFoldChange<=-0.5 & dt$Zstatistic<=-2 & dt$fdr<=0.05,]
        UPinTreat<-dt[dt$LogFoldChange>= 0.5 & dt$Zstatistic>=2& dt$fdr<=0.05,]
        print(test)
        ups<-gsub("OS\\=.*","", UPinTreat$DEFLINE)
        downs<-gsub("OS\\=.*","", UPinControl$DEFLINE)
        print("UP")
        print(ups)
        print("DOWN")
        cat(downs)
        
        upProteins[[i]]<-ups
        names(upProteins)[i]<-test
        downProteins[[i]]<-downs
        names(downProteins)[i]<-test
        
}
        
names(upProteins)
intersect(intersect(upProteins[[1]], upProteins[[2]]),upProteins[[4]])

overlaps<-intersect(upProteins[[1]], upProteins[[2]])
overlaps



#Compare the 2 results
results1<-list.files("Output/Rat/qspec_run1_coral_round/", pattern = "Results.csv")
results2<-list.files("Output/Rat/qspec_run1_coral_round_run2/", pattern = ".Results.csv")


summary<-data.frame(test=c("High.Med.All","High.Med.broadi","High","Low","Med"))
Results1<-list()
Results2<-list()
for (j in 1:2){
        if (j==1) results<-results1; path ="Output/Rat/qspec_run1_coral_round/"
        if (j==2) results<-results2; path="Output/Rat/qspec_run1_coral_round_run2/"
    
        for (i in 1: length(results)){
             dt<-read.csv(paste0(path, results[i]), stringsAsFactors = F, header=T, row.names = 1)
             UPinSite1<-dt[dt$LogFoldChange<=-0.5 & dt$Zstatistic<=-2 & dt$fdr<=0.05,]
             UPinSite2<-dt[dt$LogFoldChange>= 0.5 & dt$Zstatistic>=2& dt$fdr<=0.05,]
             
             summary[i,paste0("Up_control",j)]<-nrow(UPinSite1)
             summary[i,paste0("Up_treatment",j)]<-nrow(UPinSite2)
             if (j==1) {Results1[[i]]<-dt
                        names(Results1)[i]<-gsub(".Results.csv",'',results[i])}
             if (j==2) {Results2[[i]]<-dt
                names(Results2)[i]<-gsub(".Results.csv",'',results[i])}
             
    }
}

summary
#             test Up_control1 Up_treatment1 Up_control2 Up_treatment2
#1    High.Med.All          33             0          33             0
#2 High.Med.broadi          40             5          40             5
#3            High          52             2          52             2
#4             Low          48             3          48             3
#5             Med          57             5          57             5


#check which proteins were more abundnat in treatments
upProteins<-list()
downProteins<-list()
for (i in 1:length(Results1)){
        dt1<-Results1[[i]]
        dt2<-Results2[[i]]
        test<-names(Results1)[i]
        
        UPinControl1<-dt1[dt1$LogFoldChange<=-0.5 & dt1$Zstatistic<=-2 & dt1$fdr<=0.1,]
        UPinTreat1  <-dt1[dt1$LogFoldChange>= 0.5 & dt1$Zstatistic>=2  & dt1$fdr<=0.1,]
        
        UPinControl2<-dt2[dt2$LogFoldChange<=-0.5 & dt2$Zstatistic<=-2 & dt2$fdr<=0.1,]
        UPinTreat2  <-dt2[dt2$LogFoldChange>= 0.5 & dt2$Zstatistic>=2  & dt2$fdr<=0.1,]
        
        
        print(test)
        ups1<-gsub("OS\\=.*","", UPinTreat1$DEFLINE)
        ups2<-gsub("OS\\=.*","", UPinTreat2$DEFLINE)
        print("UP1")
        print(ups1)
        print("UP2")
        print(ups2)
        
        UPinTreat1$Zstatistic
        UPinTreat2$Zstatistic

        UPinTreat1$LogFoldChange
        UPinTreat2$LogFoldChange
        
        downs<-gsub("OS\\=.*","", UPinControl$DEFLINE)
        print("UP")
        print(ups)
        print("DOWN")
        cat(downs)
        
        upProteins[[i]]<-ups
        names(upProteins)[i]<-test
        downProteins[[i]]<-downs
        names(downProteins)[i]<-test
        
}





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


#######################################
### Rat poison exp
# look at the qspec results /attach uniprotiD/defline
results<-list.files("Output/Rat/qspec_run1_indivTest/", pattern = ".fdr")

#attach the protein annotation info
pID<-read.csv("Data/ProteinID_sorted.csv", stringsAsFactors = F) 
Results<-list()

summary<-data.frame(test=gsub(".txt_qspec_paired_fdr",'',results))
for (i in 1:length(results)){
        dt<-read.table(paste0("Output/Rat/qspec_run1_indivTest/", results[i]), stringsAsFactors = F, header=T, sep="\t")
        test<-gsub(".txt_qspec_paired_fdr", '',results[i])
        
        UPinSite1<-dt[dt$LogFoldChange<=-0.5 & dt$Zstatistic<=-2 & dt$fdr<=0.01,]
        UPinSite2<-dt[dt$LogFoldChange>= 0.5 & dt$Zstatistic>=2& dt$fdr<=0.01,]
        
        summary[i,"Up_control"]<-nrow(UPinSite1)
        summary[i,"Up_treatment"]<-nrow(UPinSite2)

        
        re1<-merge(dt,pID,by="Protein",all.x = T)
        Results[[i]]<-re1
        names(Results)[i]<-test
        write.csv(re1, paste0("Output/Rat/qspec_run1_indivTest/NonPaired_",test,".Results.csv"), row.names=T)
}

summary #30,000 iteration
#            test Up_control Up_treatment
#1 Rat_FSW.vs.HBP          7            8
#2 Rat_FSW.vs.HIP          8            1
#3 Rat_FSW.vs.LBP        133          207
#4 Rat_FSW.vs.LIP         24           20
#5 Rat_FSW.vs.MBP         36           25
#6 Rat_FSW.vs.MIP         20            2

        
summary #10,000 iteration
#            test Up_control Up_treatment
#1 Rat_FSW.vs.HBP        133           22
#2 Rat_FSW.vs.HIP         13            1
#3 Rat_FSW.vs.LBP        126           12
#4 Rat_FSW.vs.LIP        105           22
#5 Rat_FSW.vs.MBP        136           43
#6 Rat_FSW.vs.MIP         31            4
        
summary #50,000 iteration        
#            test Up_control Up_treatment
#1 Rat_FSW.vs.HBP          7            5
#2 Rat_FSW.vs.HIP          1            0
#3 Rat_FSW.vs.LBP         15            2
#4 Rat_FSW.vs.LIP          8            3
#5 Rat_FSW.vs.MBP          7            1
#6 Rat_FSW.vs.MIP        110          176

summary # run 80000 iteration with 6000 burnin for MIP (FDR=0.05)
#            test Up_control Up_treatment
#1 Rat_FSW.vs.HBP          7            5 #all fives are in lower iterations
#2 Rat_FSW.vs.HIP          1            0
#3 Rat_FSW.vs.LBP         15            2 #
#4 Rat_FSW.vs.LIP          8            3
#5 Rat_FSW.vs.MBP          7            1
#6 Rat_FSW.vs.MIP         14            4

summary #FDR=0.01
#            test Up_control Up_treatment
#1 Rat_FSW.vs.HBP          4            5
#2 Rat_FSW.vs.HIP          1            0
#3 Rat_FSW.vs.LBP         10            0
#4 Rat_FSW.vs.LIP          5            1
#5 Rat_FSW.vs.MBP          4            0
#6 Rat_FSW.vs.MIP          7            3

summary #FDR=0.1
#            test Up_control Up_treatment
#1 Rat_FSW.vs.HBP          9            5
#2 Rat_FSW.vs.HIP          1            1
#3 Rat_FSW.vs.LBP         25            7
#4 Rat_FSW.vs.LIP         13            3
#5 Rat_FSW.vs.MBP         10            4
#6 Rat_FSW.vs.MIP         20           10


upProteins<-list()
downProteins<-list()
for (i in 1:length(Results)){
        dt1<-Results[[i]]
        test<-names(Results)[i]
        
        UPinControl1<-dt1[dt1$LogFoldChange<=-0.5 & dt1$Zstatistic<=-2 & dt1$fdr<=0.01,]
        UPinTreat1  <-dt1[dt1$LogFoldChange>= 0.5 & dt1$Zstatistic>=2  & dt1$fdr<=0.01,]

        cat(paste(i,test,"\n"))
        ups<-gsub("OS\\=.*","", UPinTreat1$DEFLINE)
        #print(ups)
        #cat("\n")
        downs<-gsub("OS\\=.*","", UPinControl1$DEFLINE)
        print(downs)
        cat("\n")
        
        upProteins[[i]]<-ups
        names(upProteins)[i]<-test
        downProteins[[i]]<-downs
        names(downProteins)[i]<-test
}

#hbp<-Results[[1]]
#UPinControl1<-hbp[hbp$LogFoldChange<=-0.5 & hbp$Zstatistic<=-2 & hbp$fdr<=0.05,]
#UPinTreat1  <-hbp[hbp$LogFoldChange> 0 & hbp$fdr<=0.05,]

names(Results)
#[1] "Rat_FSW.vs.HBP" "Rat_FSW.vs.HIP" "Rat_FSW.vs.LBP" "Rat_FSW.vs.LIP"
#[5] "Rat_FSW.vs.MBP" "Rat_FSW.vs.MIP"
intersect(intersect(upProteins[[1]], upProteins[[2]]),upProteins[[4]])

overlaps<-intersect(upProteins[[1]], upProteins[[5]])
overlaps




###
#the density of null and alternative distributions.
d = read.delim("Output/Rat/qspec_run1_indivTest/Rat_FSW.vs.HBP.txt_qspec_paired_fdr", header=T, as.is=T)$Zstat
tmp = read.delim("Output/Rat/qspec_run1_indivTest/Rat_FSW.vs.HBP.txt_qspec_paired_density", header=F)
hist(d, breaks=50, xlab="Zstat", main="mydata")
ff = 500  # scaling factor: need to be adjusted in each dataset
lines(tmp$V1, tmp$V4 * ff, col=3)
lines(tmp$V1, tmp$V2 * ff, col=4)
lines(tmp$V1, tmp$V3 * ff, col=2)




### Run as non-paired samples
# look at the qspec results /attach uniprotiD/defline
results<-list.files("Output/Rat/nonpaired/", pattern = ".fdr")

#attach the protein annotation info
pID<-read.csv("Data/ProteinID_sorted.csv", stringsAsFactors = F) 
Results<-list()

summary<-data.frame(test=gsub(".txt_qspec_fdr",'',results))

for (i in 1:length(results)){
        dt<-read.table(paste0("Output/Rat/nonpaired/", results[i]), stringsAsFactors = F, header=T, sep="\t")
        test<-gsub(".txt_qspec_fdr", '',results[i])
        
        UPinSite1<-dt[dt$LogFoldChange<=-0.5 & dt$Zstatistic<=-2 & dt$fdr<=0.01,]
        UPinSite2<-dt[dt$LogFoldChange>= 0.5 & dt$Zstatistic>=2& dt$fdr<=0.01,]
        
        summary[i,"Up_control"]<-nrow(UPinSite1)
        summary[i,"Up_treatment"]<-nrow(UPinSite2)
        
        
        re1<-merge(dt,pID,by="Protein",all.x = T)
        Results[[i]]<-re1
        names(Results)[i]<-test
        write.csv(re1, paste0("Output/Rat/nonpaired/NonPaired_",test,".Results.csv"), row.names=T)
}

summary
#            test Up_control Up_treatment
#1 Rat_FSW.vs.HBP         24           25
#2 Rat_FSW.vs.HIP         21           20
#3 Rat_FSW.vs.LBP          7            8
#4 Rat_FSW.vs.LIP         32            7
#5 Rat_FSW.vs.MBP         12            8
#6 Rat_FSW.vs.MIP          8            8


upProteins<-list()
downProteins<-list()
for (i in 1:length(Results)){
        dt1<-Results[[i]]
        test<-names(Results)[i]
        
        UPinControl1<-dt1[dt1$LogFoldChange<=-0.5 & dt1$Zstatistic<=-2 & dt1$fdr<=0.01,]
        UPinTreat1  <-dt1[dt1$LogFoldChange>= 0.5 & dt1$Zstatistic>=2  & dt1$fdr<=0.01,]
        
        cat(paste(i,test,"\n"))
        ups<-gsub("OS\\=.*","", UPinTreat1$DEFLINE)
        print(ups)
        cat("\n")
        downs<-gsub("OS\\=.*","", UPinControl1$DEFLINE)
       
        
        upProteins[[i]]<-ups
        names(upProteins)[i]<-test
        downProteins[[i]]<-downs
        names(downProteins)[i]<-test
}

for (i in 1:length(Results)){
        dt1<-Results[[i]]
        test<-names(Results)[i]
        UPinControl1<-dt1[dt1$LogFoldChange<=-0.5 & dt1$Zstatistic<=-2 & dt1$fdr<=0.01,]
        UPinTreat1  <-dt1[dt1$LogFoldChange>= 0.5 & dt1$Zstatistic>=2  & dt1$fdr<=0.01,]
        cat(paste(i,test,"\n"))
        downs<-gsub("OS\\=.*","", UPinControl1$DEFLINE)
        print(downs)
        cat("\n")
}
