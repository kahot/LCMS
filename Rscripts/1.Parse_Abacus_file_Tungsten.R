##########
library(dplyr)

#
ab<-read.csv("Data/ABACUS_Plob_tungsten_output.tsv", sep = "\t",header=T)
ab<-read.csv("Data/ABACUS_Tungsten_Pdam_output.tsv", sep = "\t",header=T)


#keep only proteins that have at least 2 unique spectral counts across replicates
ab<-ab[ab$ALL_NUMPEPSUNIQ >1,] #2777

#remove contaminants
ab<-ab[!grepl("CONTAMINANT",ab$DEFLINE),]

#for P.dam
rem<-ab[!grepl("Pocillopora",ab$DEFLINE),]
ab<-ab[grepl("Pocillopora|SymbC1",ab$DEFLINE),]


# Select spectral count number columns
spec<-ab[,grepl('PROT|NUMSPECSTOT|DEFLINE', colnames(ab))]

#Rename column names by removing unnecessary strings
colnames(spec)<-gsub("_NUMSPECSTOT", "", colnames(spec))

#remove the tag for sample names
sample.name<-colnames(spec)[5]
sample.name
remove<-substr(sample.name,1,24)

n<-which(colnames(spec)=="ALL")
colnames(spec)[(n+1):ncol(spec)]<-gsub(remove, '',colnames(spec)[(n+1):ncol(spec)])

colnames(spec)[(n+1):ncol(spec)]<-substr(colnames(spec)[(n+1):ncol(spec)], 1,3)

#without pooling tech reps
write.csv(spec, "Data/Plobata_Tungsten_non-pooled.csv", row.names = F)
write.csv(spec, "Data/Pdam_Tungsten_non-pooled.csv", row.names = F)


## Pool technical replicates and use the average (if applied)

spec_pooled<-spec[,1:3]

#For P.lobata
#remove 16 (Outliers)
spec<-spec[,!(colnames(spec)=="16"|(colnames(spec)=="16B"))]
names<-colnames(spec)[(n+1):ncol(spec)]
names
names<-names[c(TRUE, FALSE)]
s<-seq(5, ncol(spec), 2)
for (i in 1:length(s)){
        samplename<-names[i]
        spec_pooled[,samplename]<-rowMeans(cbind(spec[,s[i]], spec[,(s[i]+1)]))
        
}

write.csv(spec_pooled, "Data/Plobata_tungsten.pooled.csv",row.names = F)

#For P.dam
names<-colnames(spec)[(n+1):ncol(spec)]
#remove the outliers
names<-names[!names=="33B"]

names #29B has only 1 replicate
names<-names[!names=="29B"]

spec2<-spec[,colnames(spec)%in% names]
names<-names[c(TRUE, FALSE)]

s<-seq(1, ncol(spec2), 2)
for (i in 1:length(s)){
        samplename<-names[i]
        spec_pooled[,samplename]<-rowMeans(cbind(spec2[,s[i]], spec2[,(s[i]+1)]))
}
#Attach the 29B to pooled file
spec_pooled$`29`<-spec$`29B`

write.csv(spec_pooled, "Data/Pdam_tungsten.pooled.csv",row.names = F)


#Attach Sample descriptions
samples<-read.csv("Data/2021_Oct_26_Richmond.csv")
samples<-samples[,c("File.Name", "Sample.ID")]
samples$File.Name<-gsub("2021_Oct_26_RichmondTu_",'',samples$File.Name)

spec_pooled<-read.csv("Data/Plobata_tungsten.pooled.csv",check.names = F )
spec_pooled<-read.csv("Data/Pdam_tungsten.pooled.csv",check.names = F )


samplenames<-data.frame(File.Name=colnames(spec_pooled)[4:ncol(spec_pooled)])
samplenames<-merge(samplenames, samples, by="File.Name")

samplenames
#File.Name Sample.ID
#1         11    PL30-C
#2         12   PL12-1m
#3         13 PL12-100b
#4         14 PL30-100m
#5         15 PL12-100m
#6         16 PL12-100b
#7         17  PL30-10m
#8         19 PL30-100b
#9         21   PL30-1m
#10        30    PL12-C
#11        34  PL12-10m




colnames(spec_pooled)[4:ncol(spec_pooled)]<-samplenames$Sample.ID

write.csv(spec_pooled, "Output/Plobata_tungsten.pooled.csv", row.names = F)
write.csv(spec_pooled, "Output/Pdam_tungsten.pooled.csv", row.names = F)



### Assign uniprot ID to P.lobata file
data<-read.csv("Output/Plobata_tungsten.pooled.csv", check.names = F)
colnames(data)[1]<-"Protein"
ProteinID<-read.csv("Data/ProteinID_sorted.csv", stringsAsFactors = F) 
data1<-merge(data,ProteinID[,1:2],by="Protein",all.x = T)
write.csv(data1,"Output/Tungsten/Plobata_tungsten_pooled_uniprot.csv", row.names = F)

#attach unitprot ID to Pdam file
pdam<-read.csv("Output/Pdam_tungsten.pooled.csv", check.names = F)
colnames(pdam)[1]<-"ID"
pID<-read.csv("Output/Pacuta_annotations_Uniprot.csv", stringsAsFactors = F) 
data1<-merge(pdam,pID[,c("ID","UniprotID")],by="ID",all.x = T)
write.csv(data1,"Output/Tungsten/Pdam_tungsten_pooled_uniprot.csv", row.names = F)


#FindUniprotID3<-function(data,newname){
#        data<-data
#        ProteinID<-read.csv("Data/ProteinID_sorted.csv", stringsAsFactors = F) 
#        data1<-merge(data,ProteinID[,1:3],by="Protein",all.x = T)
#        write.csv(data1,newname)
#}


#Prepare iPath file for 1 replicate per treatment

data<-read.csv("Output/Tungsten/Plobata_tungsten_pooled_uniprot.csv", )
#treat1 =Treatment 1 for comparison, treat2 =Treatment 2 for comparison
#col1 =color1, col2 #color2
treat1<-"PL30.C"
treat2<-"PL30.100m"
data[treat1]

iPATH<-function(data, treat1,treat2){
        data<-data
        ID1<-data$UniprotID[data[,treat1]!=0]
        ID1<-ID1[!is.na(ID1)]
        ID1<-paste0("UNIPROT:",ID1)
        ID2<-data$UniprotID[data[,treat2]!=0]
        ID2<-ID2[!is.na(ID2)]
        ID2<-paste0("UNIPROT:",ID2)
        write.table(ID1, paste0("Output/Tungsten/iPath.",treat1,".txt"), quote=F, sep="\t",row.names = F,col.names = F)
        write.table(ID2, paste0("Output/Tungsten/iPath.",treat2,".txt"), quote=F, sep="\t",row.names = F,col.names = F)
}

iPATH(data, "PL30.C","PL30.100m")

iPATH(data, "PL12.100m","PL30.100m")
iPATH(data, "PL12.100m","PL12.C")
iPATH(data, "PL12.1m","PL30.100b")

#For Pdam, some have 2 replicates
data<-read.csv("Output/Tungsten/Pdam_tungsten_pooled_uniprot2.csv")

#PD30.C, PD30.1m, PD30.100 have 2 replicates each 
treat<-"PD30.C"
iPATH2<-function(data,treat){
        data<-data
        df<-data[,grepl(treat, colnames(data))]
        data$Treat.sum<-apply(df, 1, sum) 
        ID1<-data$UniprotID[data$Treat.sum!=0]
        ID1<-ID1[!is.na(ID1)]
        ID1<-paste0("UNIPROT:",ID1)
        write.table(ID1, paste0("Output/Tungsten/iPath.",treat,".txt"), quote=F, sep="\t",row.names = F,col.names = F)
}
        
iPATH(data, "PD12.C","PD12.100m")
iPATH(data, "PD12.1m","PD12.10m")

iPATH2(data, "PD30.C")
iPATH2(data, "PD30.1m")
iPATH2(data, "PD30.10m")
iPATH2(data, "PD30.100m")

