# MorphoTools: a set of R functions for morphometric analysis
# 
# Plant Systematics and Evolution 301: 1115-1121, 2015
# 
# Petr Koutecký
# University of South Bohemia, Faculty of Science, Department of Botany, Branišovská 1760,
# České Budějovice, CZ-37005, Czechia
# kouta@prf.jcu.cz
#  
# Version 1.1; Mar 2016
#
# The function definitions (to be used as a source in R)
#
# ----------------------------------------------------------------------------------------
read.morphodata<-function(FILE,dec=".",popul=T){
  if (dec==".") xdata<-read.delim(FILE,header=T)
  if (dec==",") xdata<-read.delim2(FILE,header=T)
  if (popul!=T) xdata<-data.frame(xdata[,1],xdata[,2],xdata[,2],xdata[,3:ncol(xdata)])
  names(xdata)[1:3]<-c("ID","Population","Taxon")
  if (!is.factor(xdata$ID)) xdata$ID<-as.factor(xdata$ID)
  if (!is.factor(xdata$Population)) xdata$Population<-as.factor(xdata$Population)
  if (!is.factor(xdata$Taxon)) xdata$Taxon<-as.factor(xdata$Taxon)
  cat("Data strucutre using str(x):","\n")
  str(xdata)
  if (popul==T) cat("\n","List of populations: ","\n",levels(xdata$Population),"\n",fill=T)
     else cat("\n","Population column created from taxa",fill=T)
  cat("List of taxa: ","\n",levels(xdata$Taxon),fill=T)
  return(xdata)}
# ----------------------------------------------------------------------------------------
export.res<-function(RESULTS,file="clipboard",dec="."){
  write.table(RESULTS,file=file,sep="\t",quote=F,row.names=F,col.names=T,na="",dec=dec)}
# ----------------------------------------------------------------------------------------
descr.tax<-function(DATA){
  taxa<-levels(DATA$Taxon)
  noi<-seq(1,ncol(DATA)-3,1)
  res<-data.frame(factor(),factor(),numeric(),numeric(),numeric(),numeric(),numeric(),
                  numeric(),numeric(),numeric(),numeric(),numeric(),numeric())
  for (i in taxa) {
    datai<-with(DATA,DATA[Taxon==i,])
    chari<-colnames(datai)[4:ncol(datai)]
    taxi<-rep(i,ncol(datai)-3)
    lengthi<-sapply(datai[,-(1:3)],function(x) length(na.omit(x)))
    meani<-sapply(datai[,-(1:3)],mean,na.rm=T)
    sdi<-sapply(datai[,-(1:3)],sd,na.rm=T)
    quantilei<-sapply(datai[,-(1:3)],function(x) {
      quantile(x,probs=c(0.00,0.05,0.25,0.50,0.75,0.95,1.00),na.rm=T)})
    resi<-data.frame(chari,taxi,lengthi,meani,sdi,t(quantilei),noi)
    res<-rbind(res,resi)}
  names(res)<-c("Character","Taxon","N","Mean","SD","Min","5%","25%","Median","75%","95%",
                "Max","No")
  res<-with(res,res[order(No,Taxon),])
  res$Mean[which(is.nan(res$Mean))]<-NA
  res<-res[,-ncol(res)]
  return(res)}
# ----------------------------------------------------------------------------------------
descr.pop<-function(DATA){
  popul<-levels(DATA$Popul)
  noi<-seq(1,ncol(DATA)-3,1)
  res<-data.frame(factor(),factor(),numeric(),numeric(),numeric(),numeric(),numeric(),
                  numeric(),numeric(),numeric(),numeric(),numeric(),numeric())
  for (i in popul) {
    datai<-with(DATA,DATA[Population==i,])
    chari<-colnames(datai)[4:ncol(datai)]
    popi<-rep(i,ncol(datai)-3)
    lengthi<-sapply(datai[,-(1:3)],function(x) length(na.omit(x)))
    meani<-sapply(datai[,-(1:3)],mean,na.rm=T)
    sdi<-sapply(datai[,-(1:3)],sd,na.rm=T)
    quantilei<-sapply(datai[,-(1:3)],function(x) {
      quantile(x,probs=c(0.00,0.05,0.25,0.50,0.75,0.95,1.00),na.rm=T)})
    resi<-data.frame(chari,popi,lengthi,meani,sdi,t(quantilei),noi)
    res<-rbind(res,resi)}
  names(res)<-c("Character","Population","N","Mean","SD","Min","5%","25%","Median","75%",
                "95%","Max","No")
  res<-with(res,res[order(No,Population),])
  res$Mean[which(is.nan(res$Mean))]<-NA
  res<-res[,-ncol(res)]
  return(res)}
# ----------------------------------------------------------------------------------------
descr.all<-function(DATA){
  taxa<-levels(DATA$Taxon)
  res<-data.frame(factor(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),
                  numeric(),numeric(),numeric(),numeric())
  chari<-colnames(DATA)[4:ncol(DATA)]
  lengthi<-sapply(DATA[,-(1:3)],function(x) length(na.omit(x)))
  meani<-sapply(DATA[,-(1:3)],mean,na.rm=T)
  sdi<-sapply(DATA[,-(1:3)],sd,na.rm=T)
  quantilei<-sapply(DATA[,-(1:3)],function(x) quantile(x,probs=c(0.00,0.05,0.25,0.50,0.75,
                                                                 0.95,1.00),na.rm=T))
  res<-data.frame(chari,lengthi,meani,sdi,t(quantilei))
  names(res)<-c("Character","N","Mean","SD","Min","5%","25%","Median","75%","95%","Max")
  res$Mean[which(is.nan(res$Mean))]<-NA
  return(res)}
# ----------------------------------------------------------------------------------------
na.meansubst<-function(DATA){
  popul<-levels(DATA$Population)
  xdata<-DATA[1,]
  for (i in popul) {
    xp<-with(DATA,DATA[Population==i,])
    meansubst<-function(x){
      m<-mean(x,na.rm=T)
      if (is.nan(m)) m<-NA
      x[which(is.na(x))]<-m
      return(x)}
    xp[,-(1:3)]<-as.data.frame(sapply(xp[,-(1:3)],meansubst))
    xdata<-rbind(xdata,xp)}
  xdata<-xdata[-1,]
  return(xdata)}
# ----------------------------------------------------------------------------------------
popul.means<-function(DATA){
  poplist<-unique(DATA[,c(2,3)])
  xdata<-na.omit(DATA)
  xn<-aggregate(xdata[,4],by=list(Population=xdata$Population),length)
  xn<-merge(poplist,xn,all=T)
  names(xn)[3]<-"N"
  xn$N[which(is.na(xn$N))]<-0
  res<-aggregate(DATA[,-(1:3)],by=list(Population=DATA$Population),mean,na.rm=T)
  nanx<-function(x) {
    x[which(is.nan(x))]<-NA
    return(x)}
  res[,-1]<-sapply(res[,-1],nanx)
  xn<-merge(xn,res,all=T)
  return(xn)}
# ----------------------------------------------------------------------------------------
popul.otu<-function(DATA){
  xpop<-DATA[,(1:3)]
  xdata<-DATA[,-(1:3)]
  colnames(xpop)<-c("ID","Population","Taxon")
  xpop[,2:3]<-DATA[,1:2]
  xpop<-cbind(xpop,xdata)
  return(xpop)}
# ----------------------------------------------------------------------------------------
charhist<-function(DATA,taxon=NULL){
  if (is.null(taxon)) xdata<-DATA else xdata<-DATA[DATA$Taxon==taxon,]
  if (is.null(taxon)) xtax<-"all" else xtax<-taxon
  xn<-4:ncol(xdata)
    for (i in xn) {
    xh<-as.numeric(xdata[,i])
    if (is.integer(xdata[,i]) & max(xdata[,i],na.rm=T)==0) next
    if (is.integer(xdata[,i]) & max(xdata[,i],na.rm=T)==1) next
    if(sum(is.na(xh))>0) xh<-xh[-which(is.na(xh))]
    xs<-seq(min(xh),max(xh),(max(xh)-min(xh))/1000)
    hist(xh,freq=F,main="",col="grey",xlab=paste(names(xdata)[i]," (",xtax,")",sep=""))
    lines(xs,dnorm(xs,mean(xh),sqrt(var(xh))),col="red",lty=1)
    par(ask=T)}
  par(ask=F)}
# ----------------------------------------------------------------------------------------
cormat<-function (DATA,method="pearson") {
  x<-cor(DATA[,-(1:3)],use="pairwise.complete.obs",method=method)
  x<-data.frame(x)
  x<-data.frame(attr(x,"row.names"),x,row.names=NULL)
  if (method=="pearson") names(x)[1]<-"Pearson"
  if (method=="spearman") names(x)[1]<-"Spearman"
  return(x)}
# ----------------------------------------------------------------------------------------
clust<-function(DATA,method="average"){
  xdata<-DATA[,-c(2,3)]
  xdata<-data.frame(xdata,row.names=1)
  xdata<-scale(xdata)
  xdist<-dist(xdata,method="euclidean")
  if (method=="UPGMA") method<-"average"
  if (method=="ward") method<-"ward.D"
  xclust<-hclust(xdist,method=method)
  return(xclust)}
# ----------------------------------------------------------------------------------------
pca.calc<-function(DATA){
  xdata<-na.omit(DATA)
  x<-prcomp(xdata[,-(1:3)],center=T,scale.=T)
  return(x)}
# ----------------------------------------------------------------------------------------
pca.eigen<-function(RESULTS){
  x<-as.numeric(1:length(RESULTS$sdev))
  x<-sapply(x,function(x) paste("PC",x,sep=""))
  res<-sapply(RESULTS$sdev,function(x) x^2)
  names(res)<-x
  print(res,digits=3)}
# ----------------------------------------------------------------------------------------
pca.cor<-function(RESULTS,N=4){
  coeff<-data.frame(RESULTS$rotation)
  coeff<-data.frame(attr(coeff,"row.names"),coeff[,1:N],row.names=NULL)
  names(coeff)[1]<-"Character"
  sdev<-RESULTS$sdev
  sdev<-sdev[1:N]
  coeff[,-1]<-t(apply(coeff[,-1],1,function(x) x*sdev))
  return(coeff)}
# ----------------------------------------------------------------------------------------
pca.scores<-function(RESULTS,DATA,N=4){
  xdata<-na.omit(DATA)
  res<-data.frame(predict(RESULTS,xdata))
  res<-res[,1:N]
  scores<-data.frame(ID=xdata$ID,Population=xdata$Population,Taxon=xdata$Taxon,res)
  return(scores)}
# ----------------------------------------------------------------------------------------
plot.scores<-function(SCORES,axes=c(1,2),pch=NULL,col=NULL,xlab=NULL,ylab=NULL,
                      labels=FALSE,legend=FALSE,legend.pos="topright",...){
  taxlev<-levels(SCORES$Taxon)
  xtax<-as.character(SCORES$Taxon)
    replace.values<-function(x){
    for (i in 1:length(x)) xtax[which(xtax%in%taxlev[i])]<-x[i]
    return(xtax)}
  if (is.null(pch)) {pch<-16; lpch<-16} else lpch<-pch
  if (is.numeric(pch)) {if (length(pch)!=1) pch<-as.numeric(replace.values(pch))}
    else {if (length(pch)!=1) pch<-replace.values(pch)}
  if (is.null(col)) {col<-as.factor(taxlev); lcol<-col} else lcol<-col
  if (length(col)!=1) col<-replace.values(col)
  if (is.null(xlab)) xlab<-paste("axis",axes[1])
  if (is.null(ylab)) ylab<-paste("axis",axes[2])
  plot(x=SCORES[,3+axes[1]],y=SCORES[,3+axes[2]],type="n",xlab=xlab,ylab=ylab,...)
  abline(h=0,v=0,lty=2,col="grey")
  points(x=SCORES[,3+axes[1]],y=SCORES[,3+axes[2]],type="p",pch=pch,col=col,...)
  if (labels==TRUE) text(x=SCORES[,3+axes[1]],y=SCORES[,3+axes[2]],labels=SCORES$ID,
                         cex=0.5,pos=4,offset=0.5)
  if (legend==TRUE) {if (length(legend.pos)==1) legend(legend.pos,taxlev,pch=lpch,
                                                       col=lcol,bty="o")
                     else legend(legend.pos[1],legend.pos[2],taxlev,pch=lpch,col=lcol,
                                 bty="o")}}
# ----------------------------------------------------------------------------------------
point.labels<-function(DATA,names=NULL,excl=F,axes=c(1,2),pos=4,offset=0.5,...) {
  if (colnames(DATA)[1]=="ID") {
    if (is.null(names)) xdata<-DATA
    else {if (excl==FALSE) xdata<-DATA[match(names,DATA$ID),]
          else xdata<-DATA[-match(names,DATA$ID),]}
    text(x=xdata[,3+axes[1]],y=xdata[,3+axes[2]],labels=xdata$ID,pos=pos,offset=offset,...)
    }
  else {
    if (is.null(names)) xdata<-DATA
    else {if (excl==FALSE) xdata<-DATA[match(names,DATA$Character),]
          else xdata<-DATA[-match(names,DATA$Character),]}
    text(x=xdata[,1+axes[1]],y=xdata[,1+axes[2]],labels=xdata$Character,pos=pos,
         offset=offset,...)}}
# ----------------------------------------------------------------------------------------
plot.characters<-function(LOADINGS,axes=c(1,2),col="red",xlab=NULL,ylab=NULL,length=0.1,
                          xlim=c(-1,1),ylim=c(-1,1),labels=TRUE,...) {
  if (is.null(xlab)) xlab<-paste("axis",axes[1])
  if (is.null(ylab)) ylab<-paste("axis",axes[2])
  plot(x=LOADINGS[,1+axes[1]],y=LOADINGS[,1+axes[2]],type="n",xlab=xlab,ylab=ylab,
       xlim=xlim,ylim=ylim,...)
  abline(h=0,v=0,lty=2,col="grey")
  arrows(0,0,LOADINGS[,1+axes[1]],LOADINGS[,1+axes[2]],col=col,length=length,...)
  if (labels==T) text(x=LOADINGS[,1+axes[1]],y=LOADINGS[,1+axes[2]],
                      labels=LOADINGS$Character,cex=0.5,pos=4,offset=0.5)}
# ----------------------------------------------------------------------------------------
plot.biplot<-function(SCORES,LOADINGS,scale.fact=0,axes=c(1,2),pch=NULL,col.points=NULL,
                      col.arrows="red",length=0.1,xlab=NULL,ylab=NULL,labels=FALSE,
                      legend=FALSE,legend.pos="topright",...){
  if (scale.fact>0) SCORES[,-(1:3)]<-SCORES[,-(1:3)]/scale.fact
  if (scale.fact<0) LOADINGS[,-1]<-LOADINGS[,-1]*-scale.fact
  taxlev<-levels(SCORES$Taxon)
  xtax<-as.character(SCORES$Taxon)
  replace.values<-function(x){
    for (i in 1:length(x)) xtax[which(xtax%in%taxlev[i])]<-x[i]
    return(xtax)}
  if (is.null(pch)) {pch<-16; lpch<-16} else lpch<-pch
  if (is.numeric(pch)) {if (length(pch)!=1) pch<-as.numeric(replace.values(pch))}
    else {if (length(pch)!=1) pch<-replace.values(pch)}
  if (is.null(col.points)) {col.points<-as.factor(taxlev); lcol<-col.points}
    else lcol<-col.points
  if (length(col.points)!=1) col.points<-replace.values(col.points)
  if (is.null(xlab)) xlab<-paste("axis",axes[1])
  if (is.null(ylab)) ylab<-paste("axis",axes[2])
  plot(x=SCORES[,3+axes[1]],y=SCORES[,3+axes[2]],type="n",xlab=xlab,ylab=ylab,...)
  abline(h=0,v=0,lty=2,col="grey")
  points(x=SCORES[,3+axes[1]],y=SCORES[,3+axes[2]],type="p",pch=pch,col=col.points,...)
  if (labels==T) text(x=SCORES[,3+axes[1]],y=SCORES[,3+axes[2]],labels=SCORES$ID,
                      cex=0.5,pos=4,offset=0.5)
  arrows(0,0,LOADINGS[,1+axes[1]],LOADINGS[,1+axes[2]],col=col.arrows,length=length,...)
  if (labels==T) text(x=LOADINGS[,1+axes[1]],y=LOADINGS[,1+axes[2]],
                      labels=LOADINGS$Character,col=col.arrows,cex=0.5,pos=4,offset=0.5)
  if (legend==TRUE) {if (length(legend.pos)==1) legend(legend.pos,taxlev,pch=lpch,
                                                       col=lcol,bty="o")
    else legend(legend.pos[1],legend.pos[2],taxlev,pch=lpch,col=lcol,bty="o")}}
# ----------------------------------------------------------------------------------------
discr.calc<-function(DATA){
  require(vegan)
  xdata<-na.omit(DATA)
  xdata<-droplevels(xdata)
  xclass<-data.frame(model.matrix(~xdata$Taxon-1),row.names=xdata$ID)
  names(xclass)<-levels(xdata$Taxon)
  xchar<-data.frame(xdata[,-(1:3)],row.names=xdata$ID)
  discr.data<<-list("class"=xclass,"char"=xchar)
  res<-cca(discr.data$class~.,data=discr.data$char)
  return(res)}
# ----------------------------------------------------------------------------------------
discr.sum<-function(RESULTS,perm=500){
  require(vegan)
  x<-summary(RESULTS,scaling=-2)
  eig<-eigenvals(RESULTS,constrained=T)
  eig<-sapply(eig,function(x) x/(1-x))
  evar<-x$concont
  ccoeff<-spenvcor(RESULTS)
  test<-anova(RESULTS,step=perm)
  if (length(eig)>1) atest<-anova(RESULTS,step=perm,by="axis") else atest<-NULL
  res<-list(eig,evar,ccoeff,test,atest)
  names(res)<-c("eigenvalues","explained.variation","canonical.correlation.coefficients",
                "test","test.of.axes")
  return(res)}
# ----------------------------------------------------------------------------------------
discr.coef<-function(RESULTS){
  require(vegan)
  coeff<-data.frame(coef(RESULTS))
  xmean<-sapply(discr.data$char,mean)
  coeff<-data.frame(attr(coeff,"row.names"),xmean,coeff,row.names=NULL)
  colnames(coeff)[c(1,2)]<-c("Character","Mean")
  return(coeff)}
# ----------------------------------------------------------------------------------------
discr.taxa<-function(RESULTS){
  require(vegan)
  taxa<-data.frame(scores(RESULTS,display="species",scaling=-2,choices=1:RESULTS$tot.chi))
  taxa<-data.frame(attr(taxa,"row.names"),taxa,row.names=NULL)
  colnames(taxa)[1]<-"Taxon"
  return(taxa)}
# ----------------------------------------------------------------------------------------
discr.scores<-function(RESULTS,NEWDATA){
  require(vegan)
  ndata<-na.omit(NEWDATA)
  res<-predict(RESULTS,ndata,type="lc",scaling=-2)
  res<-data.frame(ID=ndata$ID,Population=ndata$Population,Taxon=ndata$Taxon,res)
  return(res)}
# ----------------------------------------------------------------------------------------
discr.bip<-function(RESULTS){
  require(vegan)
  eig<-eigenvals(RESULTS,constrained=T)
  xsum<-summary(RESULTS,axes=length(eig))
  res<-xsum$biplot
  res<-t(apply(res,1,function(x) x/(sqrt(1/(1-eig)))))
  f<-envfit(RESULTS,discr.data$char,choices=1:length(eig),permutations=0)
  f<-f$vectors[[2]]
  f<-sqrt(1/(1-f))  
  res<-res*f
  if (length(eig)>1) res<-data.frame("Character"=rownames(res),res,row.names=NULL)
  else {res<-t(res)
        res<-data.frame(res)
        res<-data.frame("Character"=row.names(res),"CCA1"=res$res,row.names=NULL)}
  return(res)}
# ----------------------------------------------------------------------------------------
discr.test<-function(RESULTS,perm=500){
  require(vegan)
  cca0<-cca(discr.data$class~1,data=discr.data$char)
  mtest<-add1(cca0,scope=formula(RESULTS),test="permutation",pstep=perm,perm.max=perm)
  mtest<-data.frame("Character"=attr(mtest,"row.names"),mtest,row.names=NULL)
  utest<-anova(RESULTS,step=perm,perm.max=perm,by="margin")
  utest<-data.frame("Character"=attr(utest,"row.names"),utest,row.names=NULL)
  xtest<-list("single.characters"=mtest,"unique.contributions"=utest)
  return(xtest)}
# ----------------------------------------------------------------------------------------
discr.step<-function(RESULTS,perm=500,p=0.05,dir="forward"){
  require(vegan)
  cca0<-cca(discr.data$class~1,data=discr.data$char)
  res<-ordistep(cca0,scope=formula(RESULTS),direction=dir,pstep=perm,perm.max=perm,
                steps=100,Pin=p)
  return(res)}
# ----------------------------------------------------------------------------------------
discr.2hist<-function(SCORES,breaks=NULL,labels=NULL,col=NULL,border=NULL,tcl=-0.5,
                      tcl.min=NULL,ylim=NULL,...){
  par(mfrow=c(1,2))
  if (is.null(breaks)) {xhist<-hist(SCORES$CCA1,plot=F); breaks<-xhist$breaks}
  taxa<-levels(SCORES$Taxon)
  if (mean(SCORES[SCORES$Taxon==taxa[1],4])>mean(SCORES[SCORES$Taxon==taxa[2],4])){
     taxa<-c(taxa[2],taxa[1])}
  hist1<-hist(SCORES[SCORES$Taxon==taxa[1],4],breaks=breaks,plot=F)
  hist2<-hist(SCORES[SCORES$Taxon==taxa[2],4],breaks=breaks,plot=F)
  if (is.null(ylim)) {ymax<-max(max(hist1$counts),max(hist2$counts))
                      ylim<-c(0,round(ymax+4,digits=-1))}
  if (is.null(col)) {col1<-NULL;col2<-NULL}
  if (length(col)==1) {col1<-col;col2<-col}
  if (length(col)==2) {col1<-col[1];col2<-col[2]}
  if (is.null(border)) {border1<-NULL;border2<-NULL}
  if (length(border)==1) {border1<-border;border2<-border}
  if (length(border)==2) {border1<-border[1];border2<-border[2]}
  if (is.null(tcl.min)) tcl.min<-tcl
  plot(hist1,col=col1,border=border1,freq=T,axes=F,ylim=ylim,main="",
       xlab=paste("canonical score (",taxa[1],")",sep=""),ylab="count",...)
  if (is.null(labels)) axis(1,at=breaks,labels=T,tcl=tcl,...)
  else {axis(1,at=breaks,labels=F,tcl=tcl.min,...)
        axis(1,at=labels,labels=as.character(labels),tcl=tcl,...)}
  axis(2,tcl=tcl,...)
  plot(hist2,col=col2,border=border2,freq=T,axes=F,ylim=ylim,main="",
       xlab=paste("canonical score (",taxa[2],")",sep=""),ylab="count",...)
  if (is.null(labels)) axis(1,at=breaks,labels=T,tcl=tcl,...)
  else {axis(1,at=breaks,labels=F,tcl=tcl.min,...)
        axis(1,at=labels,labels=as.character(labels),tcl=tcl,...)}
  axis(2,tcl=tcl,...)
  par(mfrow=c(1,1))}
# ----------------------------------------------------------------------------------------
discr.hist<-function(SCORES,breaks=NULL,labels=NULL,col=c("red","blue",0.5),tcl=-0.5,tcl.min=NULL,ylim=NULL,
                     legend=FALSE,legend.pos="topright",...){
  taxa<-levels(SCORES$Taxon)
  col2<-adjustcolor(col[2],alpha.f=col[3])
  if (is.null(breaks)) {xhist<-hist(SCORES$CCA1,plot=F); breaks<-xhist$breaks}
  hist1<-hist(SCORES[SCORES$Taxon==taxa[1],4],breaks=breaks,plot=F)
  hist2<-hist(SCORES[SCORES$Taxon==taxa[2],4],breaks=breaks,plot=F)
  if (is.null(ylim)) {ymax<-max(max(hist1$counts),max(hist2$counts))
                      ylim<-c(0,round(ymax+4,digits=-1))}
  plot(hist1,main="",xlab="canonical score",ylab="count",col=col[1],ylim=ylim,axes=F,...)
  plot(hist2,main="",xlab="canonical score",ylab="count",col=col2,ylim=ylim,
       axes=F,add=T,...)
  if (is.null(labels)) axis(1,at=breaks,labels=T,tcl=tcl,...)
  else {axis(1,at=breaks,labels=F,tcl=tcl.min,...)
        axis(1,at=labels,labels=as.character(labels),tcl=tcl,...)}
  axis(2,tcl=tcl,...)
  if (legend==TRUE) {if (length(legend.pos)==1) {
                         legend(legend.pos,levels(indivdis2.scores$Taxon),pch=22,
                                pt.bg=c(col[1],col2),pt.cex=2,bty="o")}
                    else legend(legend.pos[1],legend.pos[2],levels(indivdis2.scores$Taxon),
                                pch=22,pt.bg=c(col[1],col2),pt.cex=2,bty="o")}}
# ----------------------------------------------------------------------------------------
classif.da<-function(DATA,crossval="indiv"){
  require(MASS)
  xdata<-na.omit(DATA)
  xdata<-droplevels(xdata)
  ntax<-length(levels(xdata$Taxon))
  char<-colnames(xdata)[-c(1:3)]
  if (crossval=="ID" | crossval=="ind") crossval<-"indiv"
  if (crossval!="indiv" & crossval!="pop") stop("invalid crossvalidation unit")
  if (crossval=="indiv"){
    lda.res<-lda(as.formula(paste("Taxon ~ ",paste(char,collapse="+"))),data=xdata,
                 prior=rep(1/ntax,ntax),CV=TRUE)
    res<-data.frame(ID=xdata$ID,Population=xdata$Population,Taxon=xdata$Taxon,
                    Classif=lda.res$class,lda.res$posterior)}
  else {pop<-levels(xdata$Population)
    x<-data.frame(replicate(4,factor(),simplify=F))
    y<-data.frame(replicate(ntax,numeric(),simplify=F))
    res<-merge(x,y)
    names(res)<-c("ID","Population","Taxon","Classif",levels(xdata$Taxon))
    for (i in pop) {
      samp<-with(xdata,xdata[Population==i,])
      train<-with(xdata,xdata[(Population!=i),])
      lda.train<-lda(as.formula(paste("Taxon ~ ",paste(char,collapse="+"))),data=train,
                     prior=rep(1/ntax,ntax))
      lda.samp<-predict(lda.train,samp)
      resi<-data.frame(ID=samp$ID,Population=samp$Population,Taxon=samp$Taxon,
                       Classif=lda.samp$class,lda.samp$posterior)
      res<-rbind(res,resi)}}
  res$Correct<-as.character(res$Taxon)==as.character(res$Classif)
  return(res)}
# ----------------------------------------------------------------------------------------
classif.samp<-function(SAMPLE,TRAINING){
  require(MASS)
  samp<-na.omit(SAMPLE)
  samp<-droplevels(samp)  
  train<-na.omit(TRAINING)
  train<-droplevels(train)
  ntax<-length(levels(train$Taxon))
  char<-colnames(train)[-c(1:3)]
  lda.train<-lda(as.formula(paste("Taxon ~ ",paste(char,collapse="+"))),data=train,
                 prior=rep(1/ntax,ntax))
  lda.samp<-predict(lda.train,samp)
  res<-data.frame(ID=samp$ID,Population=samp$Population,Taxon=samp$Taxon,
                  Classif=lda.samp$class,lda.samp$posterior)
  res$Correct<-as.character(res$Taxon)==as.character(res$Classif)
  return(res)}
# ----------------------------------------------------------------------------------------
classif.matrix<-function(DATA){
  classif<-table(DATA[,3:4])
  classif<-data.frame(unclass(classif))
  classif<-data.frame(Taxon=attr(classif,"row.names"),classif,row.names=NULL)
  classif$N<-rowSums(classif[2:ncol(classif)])
  ncor<-aggregate(Correct~Taxon,data=DATA,sum)
  classif<-merge(classif,ncor)
  last<-classif[1,]
  last[1]<-"Total"
  last[2:length(last)]<-colSums(classif[2:ncol(classif)])
  classif<-rbind(classif,last)
  names(classif)[ncol(classif)]<-"percent.correct"
  classif$percent.correct<-with(classif,(percent.correct/N)*100)
  return(classif)}
# ----------------------------------------------------------------------------------------
classif.pmatrix<-function(DATA){
  classif<-table(DATA[,c(2,4)])
  classif<-data.frame(unclass(classif))
  classif<-data.frame(Population=attr(classif,"row.names"),classif,row.names=NULL)
  tax<-unique(DATA[,c(2,3)])
  classif<-merge(tax,classif)
  classif$N<-rowSums(classif[3:ncol(classif)])
  ncor<-aggregate(Correct~Population,data=DATA,sum)
  classif<-merge(classif,ncor)
  names(classif)[ncol(classif)]<-"percent.correct"
  classif$percent.correct<-with(classif,(percent.correct/N)*100)
  return(classif)}
# ----------------------------------------------------------------------------------------
classif.barplot<-function(CLASSIF,tax.ord=NULL,poplabels=NULL,col=NULL,col.lines=NULL,
                          cex.axis=1,cex.lab=1,cex.pop=1,ylab="posterior probability",
                          tcl=-0.5,ref.line=FALSE,ref.pars=list()){
  if (is.null(tax.ord)) xdata<-CLASSIF[order(CLASSIF$Taxon,CLASSIF$Population),]
  else {xdata<-CLASSIF[order(match(CLASSIF$Taxon,tax.ord),CLASSIF$Population),]
  xdata<-xdata[,c("ID","Population","Taxon","Classif",tax.ord,"Correct")]}
  if (is.null(poplabels)) poplist<-unique(xdata$Population) else poplist<-poplabels
  nind<-as.data.frame(table(xdata$Population))
  nind<-nind[order(match(nind$Var1,poplist)),]
  nind<-nind$Freq
  linepos<-c(0,cumsum(nind))
  if (is.null(col.lines)) col.lines<-par("fg")  
  barplot(as.matrix(t(xdata[,5:(ncol(xdata)-1)])),space=0,border=NA,cex.axis=cex.axis,
          cex.lab=cex.lab,col=col,ylab=ylab,width=1,xlim=c(0,sum(nind)),tcl=tcl,
          axisnames=F)
  for (i in 1:length(linepos)) abline (v=linepos[i],col=col.lines)
  textpos<-linepos[-length(linepos)]+nind/2
  mtext(poplist,at=textpos,side=1,las=2,cex=cex.pop,line=0.2)
  if (is.null(ref.pars$y)) ref.pars$y<-0.5
  if (is.null(ref.pars$col)) ref.pars$col<-par("fg")
  if (is.null(ref.pars$lwd)) ref.pars$lwd<-par("lwd")
  if (is.null(ref.pars$lty)) ref.pars$lty<-2
  if (ref.line==TRUE) lines(x=c(0,sum(nind)),y=c(ref.pars$y,ref.pars$y),col=ref.pars$col,
                            lwd=ref.pars$lwd,lty=ref.pars$lty)}
# ----------------------------------------------------------------------------------------
knn.select<-function(DATA,crossval="indiv"){
  if (crossval=="ID" | crossval=="ind") crossval<-"indiv"
  if (crossval!="indiv" & crossval!="pop") stop("invalid crossvalidation unit")
  require(class)
  k<-as.numeric(1:30)
  xdata<-na.omit(DATA)
  xdata[,-(1:3)]<-scale(xdata[,-(1:3)])
  pop<-levels(xdata$Population)
  knn.correct.1<-function(x){
    train<-xdata[,-(1:3)]      
    classif<-xdata[,3]
    knn.samp<-knn.cv(train,classif,k=x,prob=F,use.all=T)
    res<-sum(as.character(classif)==as.character(knn.samp))
    return(res)}
  knn.correct.pop<-function(x){
    res<-numeric()
    for (i in pop) {
      samp<-with(xdata,xdata[Population==i,-(1:3)])      
      samp.t<-with(xdata,xdata[Population==i,3])      
      train<-with(xdata,xdata[(Population!=i),-(1:3)])
      classif<-with(xdata,xdata[(Population!=i),3])
      knn.samp<-knn(train,samp,classif,k=x)
      resi<-sum(as.character(samp.t)==as.character(knn.samp))
      res<-sum(res,resi)}
    return(res)}
  if (crossval=="indiv") {
    for (j in 1:10){
      kselj<-sapply(k,knn.correct.1)
      if (j==1) ksel<-kselj else ksel<-rbind(ksel,kselj)}}
  else {
    for (j in 1:10){
      kselj<-sapply(k,knn.correct.pop)
      if (j==1) ksel<-kselj else ksel<-rbind(ksel,kselj)}}
  kselmean<-apply(ksel,2,mean)
  kselmax<-apply(ksel,2,max)
  kselmin<-apply(ksel,2,min)
  plot(kselmean,type="p",pch=16,xlab="K",ylab="correct classifications",
       ylim=c(min(kselmin),max(kselmax)))
  sapply(k[-1],function(x) arrows(x,kselmin[x],x,kselmax[x],code=3,angle=90,length=0.07))
  cat("The optimal K is:",which(kselmean==max(kselmean)))}
# ----------------------------------------------------------------------------------------
knn.classif<-function(DATA,K,crossval="indiv"){
  if (crossval=="ID" | crossval=="ind") crossval<-"indiv"
  if (crossval!="indiv" & crossval!="pop") stop("invalid crossvalidation unit")
  require(class)
  xdata<-na.omit(DATA)
  xdata[,-(1:3)]<-scale(xdata[,-(1:3)])
  if (crossval=="indiv") {
    train<-xdata[,-(1:3)]      
    classif<-xdata[,3]
    knn.samp<-knn.cv(train,classif,k=K,prob=T,use.all=T)
    res<-data.frame(ID=xdata$ID,Population=xdata$Population,Taxon=xdata$Taxon,
                    Classif=knn.samp,Prob=attr(knn.samp,"prob"))}
  else {
    pop<-levels(xdata$Population)
    res<-data.frame(ID=factor(),Population=factor(),Taxon=factor(),Classif=factor(),
                    Prob=numeric())
    for (i in pop) {
      samp<-with(xdata,xdata[Population==i,-(1:3)])      
      samp.a<-with(xdata,xdata[Population==i,1:3])      
      train<-with(xdata,xdata[Population!=i,-(1:3)])
      classif<-with(xdata,xdata[Population!=i,3])
      knn.samp<-knn(train,samp,classif,k=K,prob=T,use.all=T)
      resi<-data.frame(ID=samp.a$ID,Population=samp.a$Population,Taxon=samp.a$Taxon,
                       Classif=knn.samp,Prob=attr(knn.samp,"prob"))
      res<-rbind(res,resi)}}
  res$Correct<-as.character(res$Taxon)==as.character(res$Classif)
  return(res)}
# ----------------------------------------------------------------------------------------
knn.samp<-function(SAMPLE,TRAINING,K){
  require(class)
  xsamp<-na.omit(SAMPLE)
  samp<-scale(xsamp[,-(1:3)])
  xtrain<-na.omit(TRAINING)
  train<-scale(xtrain[,-(1:3)])
  classif<-xtrain[,3]
  knn.samp<-knn(train,samp,classif,k=K,prob=T,use.all=T)
  res<-data.frame(ID=xsamp$ID,Population=xsamp$Population,Taxon=xsamp$Taxon,
                  Classif=knn.samp,Prob=attr(knn.samp,"prob"))
  res$Correct<-as.character(res$Taxon)==as.character(res$Classif)
  return(res)}
