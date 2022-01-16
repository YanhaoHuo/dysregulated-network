rm(list=ls())
setwd("F:/R.3.6")
library(org.Hs.eg.db)
library(clusterProfiler)
library(edgeR)
library(openxlsx)

exp<-read.csv("exp.csv",row.names = 1)
exp1 <- exp
tpm <- log2(exp1+1)

sampletype<-read.table("sampletype.txt")[,1]
#Basal-like HER2-enriched     Luminal A     Luminal B        Normal 
sampletype <- as.vector(t(sampletype))
names(sampletype) <- colnames(exp)
LA <- names(sampletype[which(sampletype=="Luminal A")])
LB <- names(sampletype[which(sampletype=="Luminal B")])
BL <- names(sampletype[which(sampletype=="Basal-like")])
HE <- names(sampletype[which(sampletype=="HER2-enriched")])
NM <- names(sampletype[which(sampletype=="Normal")])

expLA <- tpm[,LA]
rownames(expLA) <- rownames(tpm)
expLB <- tpm[,LB]
rownames(expLB) <- rownames(tpm)
expBL <- tpm[,BL]
rownames(expBL) <- rownames(tpm)
expHE <- tpm[,HE]
rownames(expHE) <- rownames(tpm)
expNM <- tpm[,NM]
rownames(expNM) <- rownames(tpm)

IFnetwork<-read.xlsx("IFnetwork.xlsx",sheet=2)

ntwgenes <- c(t(IFnetwork[,1]),t(IFnetwork[,2]))
ntwgenes <- ntwgenes[!(duplicated(ntwgenes))]
expgenes <- intersect(rownames(tpm),ntwgenes)
notexpgs <- setdiff(ntwgenes,expgenes)

dr <- NULL
for(i in 1:nrow(IFnetwork)){
  ifrom <- intersect(IFnetwork[i,1],notexpgs)
  ito <- intersect(IFnetwork[i,2],notexpgs)
  if(length(ifrom)!=0||length(ito)!=0){
    dr <- c(dr,i)
  }
}
In <- IFnetwork[-dr,]

caledge <- function(exp,network){
  edgefrome <- exp[as.vector(t(network[,1])),]
  edgetoe <- exp[as.vector(t(network[,2])),]
  edgev <- edgefrome-edgetoe
  rownames(edgev) <- NULL
  return(edgev)
}
eLA <- caledge(expLA,In)
eLB <- caledge(expLB,In)
eBL <- caledge(expBL,In)
eHE <- caledge(expHE,In)
eNM <- caledge(expNM,In)

library(limma)
depc <- function(evn,evc,In){
  esetm <- cbind(evn,evc)
  strain <- c(rep(1,ncol(evn)),rep(2,ncol(evc)))
  design <- model.matrix(~-1+factor(strain))
  colnames(design) <- c("Normal","Tumor")
  fit <- lmFit(esetm, design)
  contrast.matrix <- makeContrasts(Tumor-Normal, levels=design) 
  fit <- eBayes(fit)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2) 
  options(digits=2)
  top<-topTable(fit2, coef=1, n=1000*1000, adjust="BH")
  avgev <- apply(evc,1,mean)
  In <- cbind(In,avgev)
  In <- In[rownames(top),]
  evc <- evc[rownames(top),]
  top <- cbind(top,In,evc)
  return(top)
}

topLA <- depc(eNM,eLA,In)
topLB <- depc(eNM,eLB,In)
topBL <- depc(eNM,eBL,In)
topHE <- depc(eNM,eHE,In)

#dysregulated network
BL<-topBL
HE<-topHE
LB<-topLB
LA<-topLA
BL_regulate<-data.frame(Gene1=BL$Gene1,Gene2=BL$Gene2,avgev=BL$avgev,logFC=BL$logFC,BL$adj.P.Val)
HE_regulate<-data.frame(Gene1=HE$Gene1,Gene2=HE$Gene2,avgev=HE$avgev,logFC=HE$logFC,HE$adj.P.Val)
LA_regulate<-data.frame(Gene1=LA$Gene1,Gene2=LA$Gene2,avgev=LA$avgev,logFC=LA$logFC,LA$adj.P.Val)
LB_regulate<-data.frame(Gene1=LB$Gene1,Gene2=LB$Gene2,avgev=LB$avgev,logFC=LB$logFC,LB$adj.P.Val)
p<-1e-4
BL_regulate<-BL_regulate[which(BL_regulate$BL.adj.P.Val<=p&(BL_regulate$logFC>=2|BL_regulate$logFC<=-2)),]
HE_regulate<-HE_regulate[which(HE_regulate$HE.adj.P.Val<=p&(HE_regulate$logFC>=2|HE_regulate$logFC<=-2)),]
LA_regulate<-LA_regulate[which(LA_regulate$LA.adj.P.Val<=p&(LA_regulate$logFC>=2|LA_regulate$logFC<=-2)),]
LB_regulate<-LB_regulate[which(LB_regulate$LB.adj.P.Val<=p&(LB_regulate$logFC>=2|LB_regulate$logFC<=-2)),]

#Key genes
source("F:/R/DS.r")#
LA_gene<-union(LA_regulate$Gene1,LA_regulate$Gene2)
LB_gene<-union(LB_regulate$Gene1,LB_regulate$Gene2)
BL_gene<-union(BL_regulate$Gene1,BL_regulate$Gene2)
HE_gene<-union(HE_regulate$Gene1,HE_regulate$Gene2)

BL_DS<-coverDS(BL_regulate,regulateDS(BL_regulate,BL_gene))
HE_DS<-coverDS(HE_regulate,regulateDS(HE_regulate,HE_gene))
LA_DS<-coverDS(LA_regulate,regulateDS(LA_regulate,LA_gene))
LB_DS<-coverDS(LB_regulate,regulateDS(LB_regulate,LB_gene))

deg<-function(BL_DS,BL_regulate){
  d<-rownames(BL_DS)
  outdeg<-c()
  indeg<-c()
  deg<-c()
  for(i in 1:length(d)){
    w1<-which(BL_regulate$Gene1==d[i])
    w2<-which(BL_regulate$Gene2==d[i])
    outdeg<-c(outdeg,length(w1))
    indeg<-c(indeg,length(w2))
    deg<-c(deg,length(w1)+length(w2))
  }
  deg<-data.frame(outdeg=outdeg,indeg=indeg,deg=deg)
  rownames(deg)<-d
  return(deg)
}
BL_DS_deg<-deg(BL_DS,BL_regulate)
HE_DS_deg<-deg(HE_DS,HE_regulate)
LA_DS_deg<-deg(LA_DS,LA_regulate)
LB_DS_deg<-deg(LB_DS,LB_regulate)

BL_CS<-rownames(BL_DS)[which(BL_DS$percentage<=0.6)]
HE_CS<-rownames(HE_DS)[which(HE_DS$percentage<=0.6)]
LA_CS<-rownames(LA_DS)[which(LA_DS$percentage<=0.6)]
LB_CS<-rownames(LB_DS)[which(LB_DS$percentage<=0.6)]


