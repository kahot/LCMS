FindUniprotID<-function(data,newname){
        data<-read.csv(data,stringsAsFactors = F)
        ProteinID<-read.csv("Data/ProteinID_sorted.csv", stringsAsFactors = F) 
        data1<-merge(data,ProteinID[,1:2],by="Protein",all.x = T)
        write.csv(data1,newname)
        
        
}