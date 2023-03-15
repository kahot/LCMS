library(reshape2)

seq<-read.table("~/programs/P.acuta/p_acuta.txt", sep="\t")
colnames(seq)<-c("ID","Seq")
anno<-read.csv("~/programs/P.acuta/Annotated/P.acuta_annotations_comb.txt", sep="\t")
colnames(anno)<-c("ID","Annotation")

Prot<-merge(seq,anno,by="ID",all.x=T)
Prot<-Prot[,c(1,3,2)]


nrow(Prot[is.na(Prot$Annotation),])#237 unannotated
Prot$Annotation[is.na(Prot$Annotation)]<-"no_annotation"
write.csv(Prot, "~/programs/P.acuta/p_acuta_proteins_annotated.csv", row.names = F)


#Merge KOALA and BLAST annotation info
koala<-read.csv("~/programs/P.acuta/split for KOALA/Done/KOALA_annotated.txt", sep="\t")
colnames(koala)<-c("ID","KO","Definition","Score","Second.best.KO","Score.second")

Prot<-read.csv("~/programs/P.acuta/p_acuta_proteins_annotated.csv")

Annotations<-merge(Prot[,1:2],koala, by="ID",all.x=T)
write.csv(Annotations,"~/programs/P.acuta/P_acuta_proteins_annotated.KO.csv" )

length(Annotations$KO[is.na(Annotations$KO)|Annotations$KO==""])/nrow(Annotations) #63.7% have no KO annotations.

#KO id list/missing list
db<-read.csv("~/programs/P.acuta/P_acuta_proteins_annotated.KO.csv", row.names = 1)
db$KOid<-db$KO
db$KOid[db$KOid==""]<-NA
db$KOid<-apply(db[c("KOid","Second.best.KO")], 1, function(x) {if (is.na(x["KOid"])) x["KOid"]<-x['Second.best.KO']
else x["KOid"]=x["KOid"]})
#Save the merge KOid file
write.csv(db,"Output/Pacuta_KO.csv")

#save the missing protiens
missing<-db[is.na(db$KOid),]
missing$Gene<-gsub("XP_\\d+.\\d\\s",'',missing$Annotation)                                                                        
missing$Gene<-gsub("NP_\\d+.\\d\\s",'',missing$Gene)                                                                        
missing$Gene<-gsub('\\[.*','',missing$Gene)                                                                        
#write.csv(missing,"Output/Pacuta_missing_KO.csv")



# Find uniprotID for annotated proteins:

db<-read.csv("Output/Pacuta_KO.csv", row.names = 1)

db$RefSeq_Name<-gsub("(XP_\\d+)(.\\d.*)",'\\1',db$Annotation)    
db$RefSeq_Name<-gsub("(NP_\\d+)(.\\d.*)",'\\1',db$RefSeq_Name)    
db$RefSeq_Name<-gsub("(WP_\\d+)(.\\d.*)",'\\1',db$RefSeq_Name)    

#write.csv(db, "Output/Pacuta.refseqnames.csv")

refseq<-unique(db$RefSeq_Name)

write.csv(db, "Output/Pacuta.refseq,unique.csv")

#Use the refseq id to find UniprotID at https://biodbnet.abcc.ncifcrf.gov/db/db2db.php#biodb
#Select Input as "RefSeq Protein Accession", Output as "Uniprot Accession"

db$Annotation<-gsub("XP_\\d+.\\d\\s",'',db$Annotation)                                                                        
db$Annotation<-gsub("NP_\\d+.\\d\\s",'',db$Annotation) 
db$Annotation<-gsub("XP_\\d+.\\d\\s",'',db$Annotation) 
db$Annotation<-gsub('\\[.*','',db$Annotation)                                                                        

uniprot<-read.csv("Data/Pacuta_Uniprot.csv")
colnames(uniprot)[1:2]<-c("RefSeq_Name","UniprotID")
uniprot<-as.data.frame(apply(uniprot,2,function(x) gsub('\\s+', '',x))) 

uniprot$UniprotID<-apply(uniprot, 1, function(x) {if (nchar(x["UniprotID"])>6) x["UniprotID"]<-x["X"]
                                                  else x["UniprotID"]<-x["UniprotID"]})

uniprot$UniprotID<-apply(uniprot, 1, function(x) {if (nchar(x["UniprotID"])>6) x["UniprotID"]<-x["X.1"]
                                                else x["UniprotID"]<-x["UniprotID"]})

write.csv(uniprot[,1:2], "Data/Pacuta_uniprot_updated.csv", row.names = F)
 
db2<-merge(db[,c("ID","Annotation","Definition","KOid","RefSeq_Name")],uniprot[,c("RefSeq_Name","UniprotID")], by="RefSeq_Name")

write.csv(db2,"Output/Pacuta_annotations_Uniprot.csv", row.names = F)
