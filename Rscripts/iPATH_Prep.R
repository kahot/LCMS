# for files with 3 replicates per treatment

iPATH_prep<-function(data, location, col1,col2,w){
        data<-read.csv(data,stringsAsFactors = F)
        if(colnames(data)[1]!="Protein") data<-data[,-1]
        data$Site1.sum<-apply(data[c(3,4,5)], 1, sum) 
        data$Site2.sum<-apply(data[c(6,7,8)], 1, sum)
        colnames(data)[9:10]<-c("LogFldChng","Zstat")
        
        ID1<-data$UniprotID[data$Site1.sum!=0]
        ID1<-ID1[!is.na(ID1)]
        ID1<-paste0("UNIPROT:",ID1)
        ID2<-data$UniprotID[data$Site2.sum!=0]
        ID2<-ID2[!is.na(ID2)]
        ID2<-paste0("UNIPROT:",ID2)
        write.table(ID1, paste0("Output/",location,"_Site1.txt"), quote=F, sep="\t",row.names = F,col.names = F)
        write.table(ID2, paste0("Output/",location,"_Site2.txt"), quote=F, sep="\t",row.names = F,col.names = F)
        
        
        UPinSite1<-data[data$LogFldChng<=-0.5 & data$Zstat<=-2,]
        UPinSite2<-data[data$LogFldChng>= 0.5 & data$Zstat>=2,]
        id1<-UPinSite1$UniprotID[!is.na(UPinSite1$UniprotID)]
        id1<-paste0("UNIPROT:",id1)
        id2<-UPinSite2$UniprotID[!is.na(UPinSite2$UniprotID)]
        id2<-paste0("UNIPROT:",id2)
        
        df1<-data.frame(uniprot=id1)
        df1$color<-col1
        df1$width<-paste0("W",w)
        df2<-data.frame(uniprot=id2)
        df2$color<-col2
        df2$width<-paste0("W",w)
        df<-rbind(df1,df2)
        
        write.table(df,paste0("Output/",location,'Site1.vs.Site2.txt'),quote=F, sep="\t",row.names = F,col.names = F)
}



# for files with 2 replicates per treatment
iPATH_prep2<-function(data, location, col1,col2,w){
        data<-read.csv(data,stringsAsFactors = F)
        if(colnames(data)[1]!="Protein") data<-data[,-1]
        data$Site1.sum<-apply(data[c(3,4)], 1, sum) 
        data$Site2.sum<-apply(data[c(5,6)], 1, sum)
        colnames(data)[7:8]<-c("LogFldChng","Zstat")
        
        ID1<-data$UniprotID[data$Site1.sum!=0]
        ID1<-ID1[!is.na(ID1)]
        ID1<-paste0("UNIPROT:",ID1)
        ID2<-data$UniprotID[data$Site2.sum!=0]
        ID2<-ID2[!is.na(ID2)]
        ID2<-paste0("UNIPROT:",ID2)
        write.table(ID1, paste0("Output/Rat/iPath.",location,"_Treat1.txt"), quote=F, sep="\t",row.names = F,col.names = F)
        write.table(ID2, paste0("Output/Rat/iPath.",location,"_Treat2.txt"), quote=F, sep="\t",row.names = F,col.names = F)
        
        
        UPinSite1<-data[data$LogFldChng<=-0.5 & data$Zstat<=-2,]
        UPinSite2<-data[data$LogFldChng>= 0.5 & data$Zstat>=2,]
        id1<-UPinSite1$UniprotID[!is.na(UPinSite1$UniprotID)]
        id1<-paste0("UNIPROT:",id1)
        id2<-UPinSite2$UniprotID[!is.na(UPinSite2$UniprotID)]
        id2<-paste0("UNIPROT:",id2)
        
        df1<-data.frame(uniprot=id1)
        df1$color<-col1
        df1$width<-paste0("W",w)
        df2<-data.frame(uniprot=id2)
        df2$color<-col2
        df2$width<-paste0("W",w)
        df<-rbind(df1,df2)
        
        write.table(df,paste0("Output/Rat/iPath_",location,'Treat1.vs.Treat2.txt'),quote=F, sep="\t",row.names = F,col.names = F)
}