library(vegan)
library(fields)
library(MASS)
library(ggplot2)
#setwd("~/Dropbox/LC-MS/NMDS")

### Palau abacus dataset ###
palau<-read.csv("~/Dropbox/LC-MS/NMDS/palau1.csv",sep=",", header = T)
#palau=apply(palau,1,as.numeric)  #This transform the matrix. Can be simply t(nnnc)

palau_nmds<-metaMDS(palau)
palau_nmds   #Stress:     5.73202e-05 
ordiplot(palau_nmds, type= "text", display="sites")  

#eliminate the first sample (too off)
palau<-palau[2:26,]
palau_nmds<-metaMDS(palau)
palau_nmds   #Stress:     0.1103268 
ordiplot(palau_nmds, type= "text", display="sites")  


#transform the data according the Brook's file
palau2<-log(palau+1)
palau_nmds1<-metaMDS(palau2, autotransform = F)
palau_nmds1   # Stress=0.04111524 
ordiplot(palau_nmds1, type= "text", display="sites")  

nikko<-c(1,2,3,4,5,6,7,8)
ng<-c(9,10,11,12,13,14,15,16)
ulong<-c(17,18,19,20,21,22,23,24,25)

ni=palau_nmds1$points[nikko,]	
nge=palau_nmds1$points[ng,]	
ul=palau_nmds1$points[ulong,]

ordiplot(palau_nmds1, type= "n", display="sites")  
points(ni, pch=21, col='red', bg='red', cex=2)
points(nge, pch=21, col='darkgreen', bg='darkgreen',cex=2)
points(ul, pch=21, col='blue', bg='blue',cex=2)


## Plot the results with ggplot
pcolors<-c("#4daf4a", "#e41a1c","#377eb8")
palau_re = as.data.frame(scores(palau_nmds1)$sites)
palau_re$Site<-c(rep("Nikko Bay", times=8),rep("Ngeremeduu Bay", times=8),rep("Ulong Island", times=9))

ggplot(palau_re, aes(x=NMDS1,y=NMDS2,color=Site)) + 
        geom_point(size=4.5) + theme_bw()+
        scale_color_manual(values=paste0(pcolors,"E6"))+
        xlim(-0.29, 0.25)+ylim(-.175,0.15)

ggsave("~/Dropbox/LC-MS/NMDS/Palau_nmds.pdf", height = 4, width = 5.7)


#calculate analysis of similarity statistics 
x<-c(rep("1", times=8),rep("2", times=8),rep("3", times=9))
ano=anosim(palau, x,permutations=1000)
ano

#anosim(x = palau, grouping = x, permutations = 1000) 
#Dissimilarity: bray 

#ANOSIM statistic R: 0.8937 
#Significance: 0.000999 

###### 

# Count the number of proteins in a GO term from COmpGo output
palbp<-read.table("~/Dropbox/LC-MS/Palau/GO/Nikko.vs.Ng/compgo_report_Nikko.vs.Ng_all_BP.txt",sep="\t", header = T)
palbp$count1<-sapply(strsplit(palbp$Protein.List..1.,","),FUN=function(x){length(x[x!="Null"])})
palbp$count2<-sapply(strsplit(palbp$Protein.List..2.,","),FUN=function(x){length(x[x!="Null"])})
write.table(palbp, "~/Dropbox/LC-MS/Palau/GO/GO_Nikkovs.Ng/compgo_report_Nikko.vs.Ng_all_BP_count.txt", sep="\t", row.names = F)

proteinCount<-function(filepath){
        df<-read.table(filepath,,sep="\t", header = T)
        df$count1<-sapply(strsplit(df$Protein.List..1.,","),FUN=function(x){length(x[x!="Null"])})
        df$count2<-sapply(strsplit(df$Protein.List..2.,","),FUN=function(x){length(x[x!="Null"])})
        write.table(df, paste0(gsub(".txt",'',filepath), "_count.txt"), sep="\t", row.names = F)
}


proteinCount("~/Dropbox/LC-MS/Palau/GO/Nikko.vs.Ng/compgo_report_Nikko.vs.Ng_all_MF.txt")
proteinCount("~/Dropbox/LC-MS/Palau/GO/Nikko.vs.Ng/compgo_report_Nikko.vs.Ng_all_CC.txt")
proteinCount("~/Dropbox/LC-MS/Palau/GO/Ng.vs.Ulong/compgo_report_Ng.vs.Ulong_all_BP.txt")
proteinCount("~/Dropbox/LC-MS/Palau/GO/Ng.vs.Ulong/compgo_report_Ng.vs.Ulong_all_MF.txt")
proteinCount("~/Dropbox/LC-MS/Palau/GO/Ng.vs.Ulong/compgo_report_Ng.vs.Ulong_all_CC.txt")


