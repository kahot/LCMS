######
# NMDS
require(vegan)
require(fields)
require(MASS)
require(fossil)
library(colorspace)
library(dplyr)
library(ggplot2)

#spec<-read.csv("Data/Brodifacoum.Alex.Oct02021_non-pooled.csv", stringsAsFactors = F,header = T)
#spec<-spec[!grepl("gi",spec$PROTID),]
#spec<-spec[!grepl("NP_",spec$PROTID),]
#spec<-spec[!grepl("P02769",spec$PROTID),]

#NMDS with NSAF
ab<-read.table("Data/AlexABACUS_Plob_poison_Oct2021_output.txt", stringsAsFactors = F,sep = "\t", header = T)

#Remove contaminant
ab<-ab[-!grepl("CONTAMINANT",ab$DEFLINE),]

#1. include all 
nsaf<-ab[,grepl('PROT|DEFLINE|ADJNSAF', colnames(ab))]

#Rename column names by removing unnecessary strings
colnames(nsaf)<-gsub("_ADJNSAF", "", colnames(nsaf))

#remove the tag for sample names
sample.name<-colnames(nsaf)[5]
remove<-substr(sample.name,1,18)

colnames(nsaf)<-gsub(remove, '',colnames(nsaf))
colnames(nsaf)<-gsub("CORAL_", '',colnames(nsaf))

df<-nsaf[-c(1:2)]
df<-t(df)
mds1<-metaMDS(df)
ordiplot(mds1, type="text", display="sites")

#unpooled data
mdsDF = as.data.frame(scores(mds1))
mdsDF$Sample<-rownames(mdsDF)
mdsDF$Sample<-substr(mdsDF$Sample, 1,3)

ggplot(mdsDF, aes(x=NMDS1,y=NMDS2,color=Sample)) + 
        geom_point(size=3.5) + theme_bw() +
        theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
        #scale_shape_manual(values=c(19,17), name="Transplant site",labels=c("Nearshore", "Offshore"))+
        xlab("Axis 1") + ylab("Axis 2")+
        theme(axis.text=element_text(size=10))
        
#Discriminant analysis with morphtools
source("Rscripts/MorphoTools1-1.R")


#Export an nsaf file in MorphoTools format
#remove the rows with all 0
df<-df[rowSums(df[,2:ncol(df)])>0,]

df<-nsaf[,c(1,3:ncol(nsaf))]
df<-t(df)
colnames(df)<-df[1,]
df<-df[-1,]

df<-data.frame(df)
df$ID<-rownames(df)
df$Population<-substr(df$ID,1,3)
df$Taxon<-df$Population
#rearrange the columns
df<-df[,c(which(colnames(df)=="ID"),which(colnames(df)=="Population"),which(colnames(df)=="Taxon"),1:3647)]

write.table(df, "Data/Alex_Ratpoison_Oct2021.NSAF.txt",row.names=F, sep ="\t" )

rat<-read.morphodata("Data/Alex_Ratpoison_Oct2021.NSAF.txt")

#
#By treatment
pop_orig<-popul.means(rat)
## populations as OTU's
populations<-popul.otu(pop_orig)


#PCA
populpca<-pca.calc(populations)
summary(populpca)
pca.eigen(populpca) # eigenvalues
plot(populpca) # scree plot

populpca.cor<-pca.cor(populpca) # character loading
populpca.scores<-pca.scores(populpca,populations) # sample scores

plot.scores(populpca.scores, legend=T,legend.pos="bottomright")

#

popd<-discr.calc(populations)



#CDA

pop.cda<-discr.calc(populations)
discr.sum(pop.cda)
#re1<-classif.da(rat,crossval="pop")

###
library(candisc)
df$Population<-factor(df$Population, levels=unique(df$Population))
vars<-colnames(df)[4:ncol(df)]
mod.rad1<- lm(cbind(vars)~Population, data=df)


#### LDA
library(MASS)
library(tidyverse)
library(caret)
data1<-df[,c(2,4:ncol(df))]
#normalize the samples  --didn't do any thing to this...

preproc.param <- data1 %>% 
        preProcess(method = c("center", "scale"))
data.transformed <- preproc.param %>% predict(data1)


#run lda on data1
model <- lda(Population~., data = data1)
plot(model)

can1<-candisc(mod1)
plot(rad1.can, scale=6)

df<-nsaf[,c(1,3:ncol(nsaf))]

y <- c(1,4,6)
d <- data.frame(y = y, x1 = c(4,-1,3), x2 = c(3,9,8), x3 = c(4,-4,-2))


#NMDS with spec count
spec<-read.csv("Data/Brodifacoum.Alex.Oct02021_non-pooled.csv", stringsAsFactors = F,header = T)


## NMDS using all proteins
#cols<-c("#F68D09","#0A41E9")        

df<-spec[-c(1:3)]
df<-t(df)
mds1<-metaMDS(df)
ordiplot(mds1, type="text", display="sites")

pn<-polanui$points[1:3,]
po<-polanui$points[4:6,]

pdf("Output/Polanui.nmds.pdf", height=6, width=6)
ordiplot(polanui, type="n", display="sites",xlim=c(-0.25,0.25))
points(pn, pch=16, col=cols[1],cex=2)
points(po, pch=16, col=cols[2],cex=2)
points(-0.23,-0.15, pch=16,col=cols[1],cex=1.5)
text(-0.20, -0.15, labels="PN",cex=1)
points(-0.23,-0.18, pch=16,col=cols[2], cex=1.5)
text(-0.20, -0.18, labels="PO",cex=1)
dev.off()




i=3
dat1<-spec %>% dplyr:: select("PROTID",grep(samples[i], names(spec)), grep(samples[i+1], names(spec)))
df<-dat1[,-1]
df<-t(df)
wahikuli<-metaMDS(df)
ordiplot(wahikuli, type="text", display="sites",xlim=c(-0.2,0.2))
wn<-wahikuli$points[1:3,]
wo<-wahikuli$points[4:6,]

pdf("Output/Wahikuli.nmds.pdf", height=6, width=6)
ordiplot(wahikuli, type="n", display="sites",xlim=c(-0.2,0.2))
points(wn, pch=16, col=cols[1],cex=2)
points(wo, pch=16, col=cols[2],cex=2)

points(-0.19,-0.12, pch=16,col=cols[1],cex=1.7)
text(-0.16, -0.12, labels="WN",cex=1)
points(-0.19,-0.14, pch=16,col=cols[2], cex=1.7)
text(-0.16, -0.14, labels="WO",cex=1)
dev.off()

i=5
dat1<-spec %>% dplyr:: select("PROTID",grep(samples[i], names(spec)), grep(samples[i+1], names(spec)))
df<-dat1[,-1]
df<-t(df)
olowalu<-metaMDS(df)
ordiplot(olowalu, type="text", display="sites")

oo<-olowalu$points[1:8,]
on<-olowalu$points[9:16,]


pdf("Output/Olowalu.nmds.pdf", height=6, width=6)
ordiplot(olowalu, type="n", display="sites",xlim=c(-0.22,0.25))
points(on, pch=16, col=cols[1],cex=2)
points(oo, pch=16, col=cols[2],cex=2)

points(-0.19,-0.13, pch=16,col=cols[1], cex=1.7)
text(-0.16, -0.13, labels="ON",cex=1)
points(-0.19,-0.15, pch=16,col=cols[2], cex=1.7)
text(-0.16, -0.15, labels="OO",cex=1)
dev.off()





##########
spec2<-spec[,-c(1:2)]
colnames(spec2)
#[1] "PN1"    "PO1"    "WN31"   "WO31"   "OO1"    "ON2"    "ON3"    "ON1"    "OO2"    "OO3"   
#[11] "OO2.1"  "WO31.1" "OO1.1"  "ON2.1"  "ON3.1"  "PN1.1"  "PO1.1"  "WN31.1" "WO31.2" "OO1.2" 
#[21] "ON2.2"  "ON3.2"  "PN1.2"  "PO1.2"  "WN31.2" "OO2.2"  "ON1.1"  "OO3.1" 

spec2<-t(spec2)

maui18<-metaMDS(spec2)
ordiplot(maui18, type="text", display="sites")

po<-maui18$points[c(2,17,24),]
wn<-maui18$points[c(3,18,25),]
wo<-maui18$points[c(4,12,19),]
pn<-maui18$points[c(1,16,23),]
on<-maui18$points[c(6:8,14,15,21,22,27),]
oo<-maui18$points[c(5,9,10,11,13,20,26,28),]

#colors<-qualitative_hcl(8,palette="Dark3")

pdf("Output/Maui2018.nmds.pdf", height=6, width=6)
xl=-0.2
xh=0.28
yl=-0.25
yh=0.15

colors<-c("#0077BB", "#EE6677", "#44BB99")
pdf("Output/Maui2018.nmds2.pdf", height=6, width=6)
ordiplot(maui18, type="n", display="sites",xlim=c(xl,xh),ylim=c(yl,yh))
points(pn, pch=16, col=colors[1],cex=1.6)
points(po, pch=17, col=colors[1],cex=1.3)
points(wn, pch=16, col=colors[2],cex=1.6)
points(wo, pch=17, col=colors[2],cex=1.3)
points(on, pch=16, col=colors[3],cex=1.6)
points(oo, pch=17, col=colors[3], cex=1.3)

points((xh-0.11),(yl+0.1), pch=16, cex=1.6)
text  ((xh-0.09),(yl+0.1), labels="Nearshore",cex=.8,adj=0)
points((xh-0.11),(yl+0.08), pch=17,cex=1.3)
text  ((xh-0.09),(yl+0.08), labels="Offshore",cex=.8, adj=0)

points((xh-0.11),(yl+0.05), pch=16,col=colors[1],cex=1.6)
text  ((xh-0.09),(yl+0.05), labels="Polanui",cex=.8, adj=0)
points((xh-0.11),(yl+0.03), pch=16,col=colors[2],cex=1.6)
text  ((xh-0.09),(yl+0.03), labels="Wahikuli",cex=.8, adj=0)
points((xh-0.11),(yl+0.01), pch=16,col=colors[3],cex=1.6)
text  ((xh-0.09),(yl+0.01), labels="Olowalu",cex=.8, adj=0)
dev.off()






######################
## NMDS using coral data only
spec1<-spec[!grepl("symb",spec$PROTID),]
spec1<-spec1[grepl("m.", spec1$PROTID),]
cols<-c("#F68D09","#0A41E9")        

i=1        
dat1<-spec1 %>% dplyr:: select("PROTID",grep(samples[i], names(spec1)), grep(samples[i+1], names(spec1)))
df<-dat1[,-1]
df<-t(df)
polanui<-metaMDS(df)
ordiplot(polanui, type="text", display="sites")

pn<-polanui$points[1:3,]
po<-polanui$points[4:6,]

pdf("Output/Polanui.nmds_coral_only.pdf", height=6, width=6)
ordiplot(polanui, type="n", display="sites",xlim=c(-0.25,0.25))
points(pn, pch=16, col=cols[1],cex=2)
points(po, pch=16, col=cols[2],cex=2)
#draw.ellipse(x=-0.06, y=-0.11, a=.12, b=.12, angle=10, border='#F8ce46', lwd=1.5)
#draw.ellipse(x=0.0, y=0.1, a=.32, b=.07, angle=158.5, border='#0433EE', lwd=1.5)
#legend(-0.4, -0.2,,legend=c("PN", "PO"), pch=16,col=cols)
points(-0.24,-0.18, pch=16,col=cols[1],cex=1.5)
text(-0.21, -0.18, labels="PN",cex=1)
points(-0.24,-0.20, pch=16,col=cols[2], cex=1.5)
text(-0.21, -0.20, labels="PO",cex=1)
dev.off()




i=3
dat1<-spec1 %>% dplyr:: select("PROTID",grep(samples[i], names(spec1)), grep(samples[i+1], names(spec1)))
df<-dat1[,-1]
df<-t(df)
wahikuli<-metaMDS(df)
ordiplot(wahikuli, type="text", display="sites",xlim=c(-0.2,0.2))

wn<-wahikuli$points[1:3,]
wo<-wahikuli$points[4:6,]

pdf("Output/Wahikuli.nmds_coral_only.pdf", height=6, width=6)
ordiplot(wahikuli, type="n", display="sites",xlim=c(-0.2,0.2))
points(wn, pch=16, col=cols[1],cex=2)
points(wo, pch=16, col=cols[2],cex=2)
points(-0.18,-0.12, pch=16,col=cols[1],cex=1.5)
text(-0.15, -0.12, labels="WN",cex=1)
points(-0.18,-0.14, pch=16,col=cols[2], cex=1.5)
text(-0.15, -0.14, labels="WO",cex=1)
dev.off()


i=5
dat1<-spec1 %>% dplyr:: select("PROTID",grep(samples[i], names(spec1)), grep(samples[i+1], names(spec1)))
df<-dat1[,-1]
df<-t(df)
olowalu<-metaMDS(df)
ordiplot(olowalu, type="text", display="sites")

oo<-olowalu$points[1:8,]
on<-olowalu$points[9:16,]


pdf("Output/Olowalu.nmds.coral.only.pdf", height=6, width=6)
ordiplot(olowalu, type="n", display="sites",xlim=c(-0.2,0.25))
points(on, pch=16, col=cols[1],cex=2)
points(oo, pch=16, col=cols[2],cex=2)

points(-0.18,-0.14, pch=16,col=cols[1],cex=1.8)
text(-0.15, -0.14, labels="ON",cex=1)
points(-0.18,-0.16, pch=16,col=cols[2],cex=1.8)
text(-0.15, -0.16, labels="OO",cex=1)
dev.off()

#WMD
df<-decostand(df,method="hellinger" )
OloWMD<-wcmdscale(vegdist(df), 2, add = TRUE,eig = TRUE)
olo_wmd1 = as.data.frame(scores(OloWMD))
ordiplot(OloWMD, type="text", display="sites")
on<-OloWMD$points[1:8,]
oo<-OloWMD$points[9:16,]

pdf("Output/Olowalu.WMS.coral.only.pdf", height=6, width=6)
xl=-0.1
xh=0.15
yl=-0.07
yh=0.13

ordiplot(OloWMD, type="n", display="sites", xlim=c(xl,xh),ylim=c(yl,yh))
points(on, pch=16, col=cols[1],cex=2)
points(oo, pch=16, col=cols[2],cex=2)
points((xh-0.05), (yh), pch=16,col=cols[1],cex=1.8)
text((xh-0.03), (yh), labels="ON",cex=1)
points((xh-0.05),(yh-0.01), pch=16,col=cols[2],cex=1.8)
text((xh-0.03), (yh-0.01), labels="OO",cex=1)
dev.off()



##########
spec2<-spec1[,-c(1:2)]
spec2H<-decostand(spec2,method="hellinger" )

#1
spec2<-t(spec2)
maui18<-metaMDS(spec2)


#2
spec2<-t(spec2H)
maui18<-metaMDS(spec2)
ordiplot(maui18, type="text", display="sites")

po<-maui18$points[c(2,17,24),]
wn<-maui18$points[c(3,18,25),]
wo<-maui18$points[c(4,12,19),]
pn<-maui18$points[c(1,16,23),]
on<-maui18$points[c(6:8,14,15,21,22,27),]
oo<-maui18$points[c(5,9,10,11,13,20,26,28),]


xl=-0.25
xh=0.3
yl=-0.2
yh=0.25

colors<-c("#0077BB", "#EE6677", "#44BB99")
pdf("Output/Maui2018.nmds.coral.only.pdf", height=6, width=6)
ordiplot(maui18, type="n", display="sites",xlim=c(xl,xh),ylim=c(yl,yh))
points(pn, pch=16, col=colors[1],cex=1.6)
points(po, pch=17, col=colors[1],cex=1.3)
points(wn, pch=16, col=colors[2],cex=1.6)
points(wo, pch=17, col=colors[2],cex=1.3)
points(on, pch=16, col=colors[3],cex=1.6)
points(oo, pch=17, col=colors[3], cex=1.3)


points((xl+0.01),(yh-0.01), pch=16, cex=1.6)
text  ((xl+0.03),(yh-0.01), labels="Nearshore",cex=.8,adj=0)
points((xl+0.01),(yh-0.03), pch=17,cex=1.3)
text  ((xl+0.03),(yh-0.03), labels="Offshore",cex=.8, adj=0)

points((xh-0.10),(yh-0.00), pch=16,col=colors[1],cex=1.6)
text  ((xh-0.08),(yh-0.00), labels="Polanui",cex=.8, adj=0)
points((xh-0.10),(yh-0.02), pch=16,col=colors[2],cex=1.6)
text  ((xh-0.08),(yh-0.02), labels="Wahikuli",cex=.8, adj=0)
points((xh-0.10),(yh-0.04), pch=16,col=colors[3],cex=1.6)
text  ((xh-0.08),(yh-0.04), labels="Olowalu",cex=.8, adj=0)

dev.off()

#Plot based on habitat types
pdf("Output/Maui2018.nmds.coral.only.by.habitats.pdf", height=6, width=6)

xl=-0.25
xh=0.3
yl=-0.2
yh=0.25

ordiplot(maui18, type="n", display="sites",xlim=c(xl,xh),ylim=c(yl,yh))
points(pn, pch=16, col=colors[2],cex=1.6)
points(po, pch=16, col=colors[1],cex=1.3)
points(wn, pch=17, col=colors[2],cex=1.6)
points(wo, pch=17, col=colors[1],cex=1.3)
points(on, pch=15, col=colors[2],cex=1.6)
points(oo, pch=15, col=colors[1], cex=1.3)

points((xl+0.01),(yh-0.01), pch=16, cex=1.6, col=colors[2])
text  ((xl+0.03),(yh-0.01), labels="Nearshore",cex=.8,adj=0)
points((xl+0.01),(yh-0.03), pch=16,cex=1.6,col=colors[1] )
text  ((xl+0.03),(yh-0.03), labels="Offshore",cex=.8, adj=0)

points((xh-0.10),(yh-0.00), pch=16,cex=1.6)
text  ((xh-0.08),(yh-0.00), labels="Polanui",cex=.8, adj=0)
points((xh-0.10),(yh-0.02), pch=17,cex=1.4)
text  ((xh-0.08),(yh-0.02), labels="Wahikuli",cex=.8, adj=0)
points((xh-0.10),(yh-0.04), pch=15,cex=1.4)
text  ((xh-0.08),(yh-0.04), labels="Olowalu",cex=.8, adj=0)
dev.off()


#### With pooled Olowalu data

dat<-spec1 %>% dplyr:: select("PROTID","DEFLINE",grep("ON", names(spec1)), grep("OO", names(spec1)))
colnames(dat)[3:18]<-substr(colnames(dat)[3:18], start=1, stop=3)
cnames<-unique(colnames(dat)[3:18])
dat2<-dat[,1:2]

for (j in 1:length(cnames)){
        dat2[,(2+j)]<-rowSums(dat[,which(colnames(dat)==cnames[j])])
        colnames(dat2)[(2+j)]<-cnames[j]
}

dat2<-dat2[,-2]
dat2S<-dat2
dat2S[,2:7]<-dat2S[,2:7]/3
spec1P<-spec1 %>% dplyr:: select("PROTID","DEFLINE",grep("PN", names(spec1)), grep("PO", names(spec1)),grep("WN", names(spec1)),grep("WO", names(spec1)))
specP<-merge(spec1P,dat2S,by="PROTID")



spec2P<-specP[,-c(1:2)]
spec2H<-decostand(spec2P,method="hellinger" )

#2
spec2<-t(spec2H)
maui18<-metaMDS(spec2)
ordiplot(maui18, type="text", display="sites")

rownames(spec2)
po<-maui18$points[c(4:6),]
wn<-maui18$points[c(7:9),]
wo<-maui18$points[c(10:12),]
pn<-maui18$points[c(1:3),]
on<-maui18$points[c(13:15),]
oo<-maui18$points[c(16:18),]

xl=-0.3
xh=0.3
yl=-0.2
yh=0.25

ordiplot(maui18, type="n", display="sites",xlim=c(xl,xh),ylim=c(yl,yh))
points(pn, pch=16, col=colors[2],cex=1.6)
points(po, pch=16, col=colors[1],cex=1.3)
points(wn, pch=17, col=colors[2],cex=1.6)
points(wo, pch=17, col=colors[1],cex=1.3)
points(on, pch=15, col=colors[2],cex=1.6)
points(oo, pch=15, col=colors[1], cex=1.3)

points((xl+0.3),(yh-0.01), pch=16, cex=1.6, col=colors[2])
text  ((xl+0.32),(yh-0.01), labels="Nearshore",cex=.8,adj=0)
points((xl+0.3),(yh-0.03), pch=16,cex=1.6,col=colors[1] )
text  ((xl+0.32),(yh-0.03), labels="Offshore",cex=.8, adj=0)

points((xh-0.10),(yh-0.00), pch=16,cex=1.6)
text  ((xh-0.08),(yh-0.00), labels="Polanui",cex=.8, adj=0)
points((xh-0.10),(yh-0.02), pch=17,cex=1.4)
text  ((xh-0.08),(yh-0.02), labels="Wahikuli",cex=.8, adj=0)
points((xh-0.10),(yh-0.04), pch=15,cex=1.4)
text  ((xh-0.08),(yh-0.04), labels="Olowalu",cex=.8, adj=0)
dev.off()


nmdsP = as.data.frame(scores(maui18))

samplenames<-rownames(nmdsP)
samplenames<-substr(samplenames, start=1, stop=2)
Samples<-substr(samplenames, start=1,stop=1)
for (i in 1:length(Samples)){
        if (Samples[i]=="P") Samples[i]<-"Polanui"
        if (Samples[i]=="W") Samples[i]<-"Wahikuli"
        if (Samples[i]=="O") Samples[i]<-"Olowalu"
}

locations<-substr(samplenames, start=2,stop=2)
for (i in 1:length(locations)){
        if (locations[i]=="N") locations[i]<-"Nearshore"
        if (locations[i]=="O") locations[i]<-"Offshore"
}



ggplot(nmdsP, aes(x=NMDS1,y=NMDS2,color=Samples,shape=locations)) + 
        geom_point(size=4) + theme_bw() +
        theme(legend.text = element_text(size=14), legend.title = element_text(size=16,face = "bold"))+
        labs(color="Locations",shape="")+
        theme(axis.text=element_text(size=10))
ggsave("./Output/Maui2018_nmds.pooledOlowalu.pdf", height = 8, width = 9.8)





###################
# Weighted Metric multi-dimensional scaling ####

wcmd1 = (wcmdscale(vegdist(spec2), 2, add = TRUE,eig = TRUE))
maui_wcmd1 = as.data.frame(scores(wcmd1))

samplenames<-rownames(maui_wcmd1)
samplenames<-substr(samplenames, start=1, stop=2 )
samples<-substr(samplenames, start=1,stop=1)
for (i in 1:length(samples)){
        if (samples[i]=="P") samples[i]<-"Polanui"
        if (samples[i]=="W") samples[i]<-"Wahikuli"
        if (samples[i]=="O") samples[i]<-"Olowalu"
}

locations<-substr(samplenames, start=2,stop=2)
for (i in 1:length(locations)){
        if (locations[i]=="N") locations[i]<-"Nearshore"
        if (locations[i]=="O") locations[i]<-"Offshore"
}



ggplot(maui_wcmd1, aes(x=Dim1,y=Dim2,color=samples,shape=locations)) + 
        geom_point(size=4) + theme_bw() +
        theme(legend.text = element_text(size=14), legend.title = element_text(size=16,face = "bold"))+
               labs(color="Locations",shape="")+
        xlab("Axis 1") + ylab("Axis 2")+
        theme(axis.text=element_text(size=10))
ggsave("./Output/Maui2018_WMD.pdf", height = 8, width = 9.8)



maui_wcmd1.2<-cbind(maui_wcmd1, samplenames)
ordiplot(maui_wcmd1, type="n", display="sites", xlab="Axis 1", ylab= "Axis 2")
points(maui_wcmd1.2[maui_wcmd1.2$samplenames=="PN",1:2],pch=16, col=colors[1],cex=1.6)
points(maui_wcmd1.2[maui_wcmd1.2$samplenames=="PO",1:2],pch=17, col=colors[1],cex=1.3)
points(maui_wcmd1.2[maui_wcmd1.2$samplenames=="WN",1:2],pch=16, col=colors[7],cex=1.6)
points(maui_wcmd1.2[maui_wcmd1.2$samplenames=="WO",1:2],pch=17, col=colors[7],cex=1.3)
points(maui_wcmd1.2[maui_wcmd1.2$samplenames=="ON",1:2],pch=16, col=colors[6],cex=1.6)
points(maui_wcmd1.2[maui_wcmd1.2$samplenames=="OO",1:2],pch=17, col=colors[6],cex=1.3)



##########################
##########################
### Use NSAF to do nmds

dt<-read.table("Data/Alex_ABACUS_output.tsv", stringsAsFactors = F,sep = "\t", header = T)
ID<-read.csv("Data/alex_sequence_file.csv", stringsAsFactors = F)
ID$No<-gsub("2019_April_12_coral_","",ID$File.Name)

#Select Adjusted NSAF value column
nsaf.dt<-dplyr::select(dt, contains("_ADJNSAF"))
nsaf.names<-gsub("X2019_APRIL_12_BACKMAN_CORAL_","", names(nsaf.dt))
nsaf.name2<-gsub("_ADJNSAF", "" , nsaf.names)

sampleNames<-data.frame(ID=rep(NA,length(nsaf.name2)))
for (i in 1:length(nsaf.name2)){
        sampleNames$ID[i]<-ID$Comment[ID$No==nsaf.name2[i]]
}

# Select coral protein only
colnames(nsaf.dt)<- sampleNames$ID
nsaf.dt<-cbind(dt[,1:2], nsaf.dt)
nsaf.dt<-nsaf.dt[!grepl("gi",nsaf.dt$PROTID),]
nsaf.dt<-nsaf.dt[!grepl("symb",nsaf.dt$PROTID),]

#
nsaf2<-nsaf.dt[,-c(1:2)]
colnames(nsaf2)
#[1] "PN1"    "PO1"    "WN31"   "WO31"   "OO1"    "ON2"    "ON3"    "ON1"    "OO2"    "OO3"   
#[11] "OO2.1"  "WO31.1" "OO1.1"  "ON2.1"  "ON3.1"  "PN1.1"  "PO1.1"  "WN31.1" "WO31.2" "OO1.2" 
#[21] "ON2.2"  "ON3.2"  "PN1.2"  "PO1.2"  "WN31.2" "OO2.2"  "ON1.1"  "OO3.1" 

nsaf2<-t(nsaf2)
nsaf3<-decostand(nsaf2,method="hellinger")
NSAFmaui<-metaMDS(nsaf3)
ordiplot(NSAFmaui, type="text", display="sites")

po<-NSAFmaui$points[c(2,17,24),]
wn<-NSAFmaui$points[c(3,18,25),]
wo<-NSAFmaui$points[c(4,12,19),]
pn<-NSAFmaui$points[c(1,16,23),]
on<-NSAFmaui$points[c(6:8,14,15,21,22,27),]
oo<-NSAFmaui$points[c(5,9,10,11,13,20,26,28),]

#colors<-qualitative_hcl(8,palette="Dark3")

xl=-0.18
xh=0.18
yl=-0.15
yh=0.15


pdf("Output/Maui2018_NSAF.nmds.pdf", height=6, width=6)
ordiplot(NSAFmaui, type="n", display="sites", xlim=c(xl,xh))
points(pn, pch=16, col=colors[1],bg= "#F8CE46", cex=1.6,lwd=2)
points(po, pch=17, col=colors[1],bg= "#0433FFB3", cex=1.3,lwd=2)
points(wn, pch=16, col=colors[2],bg= "#F8CE46", cex=1.6,lwd=2)
points(wo, pch=17, col=colors[2],bg= "#0433FFB3", cex=1.3,lwd=2)
points(on, pch=16, col=colors[3],bg= "#F8CE46", cex=1.6,lwd=2)
points(oo, pch=17, col=colors[3],bg= "#0433FFB3", cex=1.3,lwd=2)

points((xh-0.1), (yl+0.04), pch=16, col=colors[1], cex=1.6)
text  ((xh-0.08), (yl+0.04), labels="PN",cex=.8)
points((xh-0.05), (yl+0.04), pch=17,col=colors[1],cex=1.3)
text  ((xh-0.03), (yl+0.04), labels="PO",cex=.8)
points((xh-0.1), (yl+0.025), pch=16,col=colors[2],cex=1.6)
text  ((xh-0.08),(yl+0.025), labels="WN",cex=.8)
points((xh-0.05),(yl+0.025), pch=17,col=colors[2],cex=1.3)
text  ((xh-0.03),(yl+0.025), labels="WO",cex=.8)
points((xh-0.1) ,(yl+0.01), pch=16,col=colors[3],cex=1.6)
text  ((xh-0.08),(yl+0.01), labels="ON",cex=.8)
points((xh-0.05),(yl+0.01), pch=17,col=colors[3],cex=1.3)
text  ((xh-0.03),(yl+0.01), labels="OO",cex=.8)
dev.off()

        
        

###########
olo<-read.table("Output/OOvsON_pooled_corals_adj.txt",sep="\t",stringsAsFactors = F, header=T)
olo1<-olo[,c(3:8)]
olo1t<-t(olo1)
oloH<-decostand(olo1,method="hellinger")
oloHt<-t(oloH)


olo.adjusted<-olo1
olo.adjusted$S0.2<-olo.adjusted$S0.2/2*3
olo.adjusted$S1.2<-olo.adjusted$S1.2/2*3
oloH2<-decostand(olo.adjusted, method="hellinger")
olo2<-t(oloH2)
olo1<-t(oloH)

Olonmds<-metaMDS(oloHt)
ordiplot(Olonmds, type='text', display="sites")


Olonmds2<-metaMDS(olo2)
ordiplot(Olonmds2, type='text', display="sites")


oo<-Olonmds$points[1:3,]
on<-Olonmds$points[4:6,]

pdf("Output/Olowalu.pooled.coral.adj.nmds.pdf", height=6, width=6)
ordiplot(Olonmds, type="n", display="sites")
points(on, pch=16, col=cols[1],cex=2)
points(oo, pch=16, col="#0A41E9CC",cex=2)
dev.off()





# Weighted Metric multi-dimensional scaling ####

nsaf3<-decostand(nsaf2,method="hellinger")
wcmd2 = (wcmdscale(vegdist(nsaf3), 2, add = TRUE,eig = TRUE))
maui_wcmd2 = as.data.frame(scores(wcmd2))

samplenames<-rownames(maui_wcmd2)
samplenames<-substr(samplenames, start=1, stop=2 )
samples<-substr(samplenames, start=1,stop=1)
for (i in 1:length(samples)){
        if (samples[i]=="P") samples[i]<-"Polanui"
        if (samples[i]=="W") samples[i]<-"Wahikuli"
        if (samples[i]=="O") samples[i]<-"Olowalu"
}

locations<-substr(samplenames, start=2,stop=2)
for (i in 1:length(locations)){
        if (locations[i]=="N") locations[i]<-"Nearshore"
        if (locations[i]=="O") locations[i]<-"Offshore"
}

ggplot(maui_wcmd2, aes(x=Dim1,y=Dim2,color=samples,shape=locations)) + 
        geom_point(size=4) + theme_bw() +
        theme(legend.text = element_text(size=14), legend.title = element_text(size=16,face = "bold"))+
        labs(color="Locations",shape="")+
        xlab("Axis 1") + ylab("Axis 2")+
        theme(axis.text=element_text(size=10))+
        xlim(-0.2,0.18)+ylim(-0.08, 0.08)
ggsave("./Output/Maui2018_NSAF.WMD.pdf", height = 8, width = 9.8)

cols<-c("#F68D09CC","#0A41E9CC")        
ggplot(maui_wcmd2, aes(x=Dim1,y=Dim2,color=locations,shape=samples)) + 
        geom_point(size=4) + theme_bw() +
        theme(legend.text = element_text(size=14), legend.title = element_text(size=16,face = "bold"))+
        labs(color="",shape="")+
        xlab("Axis 1") + ylab("Axis 2")+
        theme(axis.text=element_text(size=10))+
        xlim(-0.2,0.18)+ylim(-0.08, 0.08)+
        scale_color_manual(values=cols)
ggsave("./Output/Maui2018_NSAF.WMD.by.habitats.pdf", height = 8, width = 9.8)


### NMDS by habitats
maui_nmds = as.data.frame(scores(NSAFmaui))
#cols<-c("#F68D09CC","#0A41E9CC")        
ggplot(maui_nmds, aes(x=NMDS1,y=NMDS2,color=locations,shape=samples)) + 
        geom_point(size=4) + theme_bw() +
        theme(legend.text = element_text(size=14), legend.title = element_text(size=16,face = "bold"))+
        labs(color="",shape="")+
        xlab("Axis 1") + ylab("Axis 2")+
        theme(axis.text=element_text(size=10))+
        xlim(-0.2,0.18)+ylim(-0.08, 0.08)+
        scale_color_manual(values=cols)
ggsave("./Output/Maui2018_NSAF.nmds.by.habitats.pdf", height = 8, width = 9.8)


