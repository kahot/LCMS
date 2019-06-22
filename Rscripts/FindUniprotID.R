FindUniprotID<-function(data,newname){
        data<-read.table(data, sep="\t", stringsAsFactors = F,  header=T)
        #data<-read.csv(data,stringsAsFactors = F)
        ProteinID<-read.csv("Data/ProteinID_sorted.csv", stringsAsFactors = F) 
        data1<-merge(data,ProteinID[,1:3],by="Protein",all.x = T)
        write.csv(data1,newname)
        
        
}

