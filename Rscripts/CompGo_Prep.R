#for 3 replicates per treatment/site
CompGo_prep<-function(data, location, col1,col2,w){
        data<-read.csv(data,stringsAsFactors = F)
        if(colnames(data)[1]!="Protein") data<-data[,-1]
        data$Site1.sum<-apply(data[c(3,4,5)], 1, sum) 
        data$Site2.sum<-apply(data[c(6,7,8)], 1, sum)
        colnames(data)[9:10]<-c("LogFldChng","Zstat")
        
        ID1<-data$Protein[data$Site1.sum!=0]
        ID1<-ID1[!is.na(ID1)]
        ID2<-data$Protein[data$Site2.sum!=0]
        ID2<-ID2[!is.na(ID2)]
        
        write.table(ID1, paste0("Output/Rat/CompGo.",location,"_Treat1.txt"), quote=F, sep="\t",row.names = F,col.names = F)
        write.table(ID2, paste0("Output/Rat/CompGo.",location,"_Treat2.txt"), quote=F, sep="\t",row.names = F,col.names = F)
        
        
        UPinSite1<-data[data$LogFldChng<=-0.5 & data$Zstat<=-2,]
        UPinSite2<-data[data$LogFldChng>= 0.5 & data$Zstat>=2,]
        id1<-UPinSite1$Protein[!is.na(UPinSite1$Protein)]
        id2<-UPinSite2$Protein[!is.na(UPinSite2$Protein)]
        
        df1<-data.frame(Protein=id1)
        df1$treat<-"Treat1"
        df2<-data.frame(Protein=id2)
        df2$treat<-"Treat2"
        df<-rbind(df1,df2)
        
        write.table(df,paste0("Output/Rat/CompGo_",location,'Treat1.vs.Treat2.txt'),quote=F, sep="\t",row.names = F,col.names = F)
}


#for 2 replicates per treatment/site
CompGo_Prep2<-function(data, location){
        data<-read.csv(data,stringsAsFactors = F)
        if(colnames(data)[1]!="Protein") data<-data[,-1]
        data$Site1.sum<-apply(data[c(3,4)], 1, sum) 
        data$Site2.sum<-apply(data[c(5,6)], 1, sum)
        colnames(data)[7:8]<-c("LogFldChng","Zstat")
        
        ID1<-data$Protein[data$Site1.sum!=0]
        ID1<-ID1[!is.na(ID1)]
        
        ID2<-data$Protein[data$Site2.sum!=0]
        ID2<-ID2[!is.na(ID2)]
        
        write.table(ID1, paste0("Output/Rat/CompGo.",location,"_Treat1.txt"), quote=F, sep="\t",row.names = F,col.names = F)
        write.table(ID2, paste0("Output/Rat/CompGo.",location,"_Treat2.txt"), quote=F, sep="\t",row.names = F,col.names = F)
        
        
        UPinSite1<-data[data$LogFldChng<=-0.5 & data$Zstat<=-2,]
        UPinSite2<-data[data$LogFldChng>= 0.5 & data$Zstat>=2,]
        id1<-UPinSite1$Protein[!is.na(UPinSite1$Protein)]
        
        id2<-UPinSite2$Protein[!is.na(UPinSite2$Protein)]
        
        
        df1<-data.frame(Protein=id1)
        df1$treat<-"Treat1"
        df2<-data.frame(Protein=id2)
        df2$treat<-"Treat2"
        df<-rbind(df1,df2)
        
        write.table(df,paste0("Output/Rat/CompGo_",location,'Treat1.vs.Treat2.txt'),quote=F, sep="\t",row.names = F,col.names = F)
}
