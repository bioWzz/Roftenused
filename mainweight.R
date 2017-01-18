rm(list=ls())
library(GenomicFeatures)
library(INSPEcT)
library(Rsamtools)
library(deSolve)
library(BiocParallel)
library(compiler)
library(sqldf)
# load("D:/labAandC.RData")
LabArep1E=GetExprssionData(LabArep1)
LabArep2E=GetExprssionData(LabArep2)
DrawingIntegerScatterPlot(LabArep1E,LabArep2E,12,1,12,"LabA.png",3,4,"Expression",c(0,1000),c(0,1000),"red")
# normalization.factors <- read.delim("C:/Users/wzz/Desktop/HarmExampleData/analysis_nobatch_FDR_q0.05/normalization-factors.txt")
LabArep1Weight=c(1.1067760,1.0939119,0.9499632,1.1067760,1.0939119,0.9499632,1.0643657,1.0575026,0.7806148,1.0643657,1.0575026,0.7806148)
LabArep2Weight=c(1.0842723,1.0900224,0.9739547,1.0842723,1.0900224,0.9739547,1.0684246,1.0530198,0.7055990,1.0684246,1.0530198,0.7055990)

LabArep1EaW=ChangeWeightsOfExpression(LabArep1E,LabArep1Weight)
LabArep2EaW=ChangeWeightsOfExpression(LabArep2E,LabArep2Weight)
DrawingIntegerScatterPlot(LabArep1EaW,LabArep2EaW,12,1,12,"LabAW.png",3,4,"Expression",c(0,1000),c(0,1000),"red")

#COMPARE EXPRESSION WIHT HARM
LabArep1E=getvaluelog2(LabArep1)
LabArep2E=getvaluelog2(LabArep2)
HLabArep1E=logCPM[,c(1,2,3,7,8,9)]
HLabArep2E=logCPM[,c(13,14,15,19,20,21)]


HLabArep1Emerge=mergeandchooseOneData(LabArep1E,HLabArep1E,2)
LabArep1Emerge=mergeandchooseOneData(LabArep1E,HLabArep1E,1)
colnames(LabArep1Emerge)=colnames(HLabArep1Emerge)
DrawingIntegerScatterPlot(LabArep1Emerge,HLabArep1Emerge,6,1,6,"LabAhaREP1.png",2,3,"Expression",c(-1,20),c(-1,20),"red")


HLabArep2Emerge=mergeandchooseOneData(LabArep2E,HLabArep2E,2)
LabArep2Emerge=mergeandchooseOneData(LabArep2E,HLabArep2E,1)
colnames(LabArep2Emerge)=colnames(HLabArep2Emerge)
DrawingIntegerScatterPlot(LabArep2Emerge,HLabArep2Emerge,6,1,6,"LabAhaREP2.png",2,3,"Expression",c(-1,20),c(-1,20),"red")


LabCrep1E=getvaluelog2(LabCrep1)
LabCrep2E=getvaluelog2(LabCrep2)
HLabCrep1E=logCPM[,c(4,5,6,10,11,12)]
HLabCrep2E=logCPM[,c(16,17,18,22,23,24)]


HLabCrep1Emerge=mergeandchooseOneData(LabCrep1E,HLabCrep1E,2)
LabCrep1Emerge=mergeandchooseOneData(LabCrep1E,HLabCrep1E,1)
colnames(LabCrep1Emerge)=colnames(HLabCrep1Emerge)
DrawingIntegerScatterPlot(LabCrep1Emerge,HLabCrep1Emerge,6,1,6,"LabChaREP1.png",2,3,"Expression",c(-1,20),c(-1,20),"red")

HLabCrep2Emerge=mergeandchooseOneData(LabCrep2E,HLabCrep2E,2)
LabCrep2Emerge=mergeandchooseOneData(LabCrep2E,HLabCrep2E,1)
colnames(LabCrep2Emerge)=colnames(HLabCrep2Emerge)
DrawingIntegerScatterPlot(LabCrep2Emerge,HLabCrep2Emerge,6,1,6,"LabChaREP2.png",2,3,"Expression",c(-1,20),c(-1,20),"red")



# subset from the alex
subset <- read.csv("D:/ivan/ivan6-cell/subset/subset.txt", sep="")
colnames(subset)="id"
load("C:/Users/wzz/Desktop/harmlogCPM.RData")
LimmaResult=ConvertENSGtoGenenameNoRownames(logCPM)
subsetOLimmaResult=merge(subset,LimmaResult,by="id")
rownames(subsetOLimmaResult)=subsetOLimmaResult[,1]
subsetOLimmaResult=subsetOLimmaResult[,2:(length(subsetOLimmaResult[1,]))]

# limma rep1 and rep2的相关性
HLabArep1E=subsetOLimmaResult[,c(1,2,3,7,8,9)]
HLabArep2E=subsetOLimmaResult[,c(13,14,15,19,20,21)]
DrawingIntegerScatterPlot(HLabArep1E,HLabArep2E,6,1,6,"LimmaLabAREP12.png",2,3,"Expression",c(-10,20),c(-10,20),"red")

HLabCrep1E=subsetOLimmaResult[,c(4,5,6,10,11,12)]
HLabCrep2E=subsetOLimmaResult[,c(16,17,18,22,23,24)]
DrawingIntegerScatterPlot(HLabCrep1E,HLabCrep2E,6,1,6,"LimmaLabCREP12.png",2,3,"Expression",c(-10,20),c(-10,20),"red")

# inspect rep1 and rep2的相关性
LabArep1E=getvaluelog2subset(LabArep1,subset)
LabArep2E=getvaluelog2subset(LabArep2,subset)
DrawingIntegerScatterPlot(LabArep1E,LabArep2E,6,1,6,"InspectLabAREP12.png",2,3,"Expression",c(-10,20),c(-10,20),"red")

LabCrep1E=getvaluelog2subset(LabCrep1,subset)
LabCrep2E=getvaluelog2subset(LabCrep2,subset)
DrawingIntegerScatterPlot(LabCrep1E,LabCrep2E,6,1,6,"InspectLabCREP12.png",2,3,"Expression",c(-10,20),c(-10,20),"red")

# inspect 和limma 之间的相关性
LabArep1Emerge=mergeandchooseOneData(LabArep1E,HLabArep1E,1)
HLabArep1Emerge=mergeandchooseOneData(LabArep1E,HLabArep1E,2)
DrawingIntegerScatterPlot(LabArep1Emerge,HLabArep1Emerge,6,1,6,"subsetLabAhaREP1.png",2,3,"Expression",c(-10,10),c(-10,10),"red")

LabArep2Emerge=mergeandchooseOneData(LabArep2E,HLabArep2E,1)
HLabArep2Emerge=mergeandchooseOneData(LabArep2E,HLabArep2E,2)
DrawingIntegerScatterPlot(LabArep2Emerge,HLabArep2Emerge,6,1,6,"subsetLabAhaREP2.png",2,3,"Expression",c(-10,10),c(-10,10),"red")

LabCrep1Emerge=mergeandchooseOneData(LabCrep1E,HLabCrep1E,1)
HLabCrep1Emerge=mergeandchooseOneData(LabCrep1E,HLabCrep1E,2)
DrawingIntegerScatterPlot(LabCrep1Emerge,HLabCrep1Emerge,6,1,6,"subsetLabChaREP1.png",2,3,"Expression",c(-10,10),c(-10,10),"red")

LabCrep2Emerge=mergeandchooseOneData(LabCrep2E,HLabCrep2E,1)
HLabCrep2Emerge=mergeandchooseOneData(LabCrep2E,HLabCrep2E,2)
DrawingIntegerScatterPlot(LabCrep2Emerge,HLabCrep2Emerge,6,1,6,"subsetLabChaREP2.png",2,3,"Expression",c(-10,10),c(-10,10),"red")

#limma，inpect得到的结果与所给集合的交集
allSET=as.character(subset[,1]) 
limmaSET=rownames(subsetOLimmaResult)
InpsectSET=rownames(subsetResult)

notinlimma=setdiff(allSET,limmaSET)
notininpsect=setdiff(allSET,InpsectSET)
notininpsect=setdiff(limmaSET,InpsectSET)

# du ru biao da shuju
# du wen jian 
for(name in c("A","C")){
    for(hour in c("4","8")){
     assign(paste("lab",name,"_",hour,sep=""),read.delim(paste("D:/ivan/ivan6-cell/subset/express/Lab",name,"_",hour,"_0.txt",sep=""))[,c(1,6)]) 
     z=get(paste("lab",name,"_",hour,sep=""));colnames(z)=c("id",paste("lab",name,"_",hour,sep=""));assign(paste("lab",name,"_",hour,sep=""),z)

     assign(paste("Unl",name,"_",hour,sep=""),read.delim(paste("D:/ivan/ivan6-cell/subset/express/Unl",name,"_",hour,"_0.txt",sep=""))[,c(1,6)]) 
     p=get(paste("Unl",name,"_",hour,sep=""));colnames(p)=c("id",paste("Unl",name,"_",hour,sep=""));assign(paste("Unl",name,"_",hour,sep=""),p)
     }
}
# he bing wen jian
z=as.data.frame(labA_4[,1]);colnames(z)=c("id")
for(name in c("A","C")){
 for(type in c("lab","Unl")){
   for(hour in c("4","8"))
   {
    z=merge(z,get(paste(type,name,"_",hour,sep="")),by="id")
   }
 }
}
# find subset gene exprssion of limma
subsetgeneExpression=merge(subset,ConvertENSGtoGenenameNoRownames(ConverttoDataWithRownames(z)),by="id")
write.table(subsetgeneExpression, file = "D:/ivan/ivan6-cell/figureResult/alex/limma_subset_geneexpression/limma_subset_geneexpression.txt", append = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

#find subset gene expression of inspect
LabArep1E=GetExprssionData(LabArep1)[,c(1,2,3,7,8,9)]
colnames(LabArep1E)=c("LabA_0","LabA_4","LabA_8","UnlA_0","UnlA_4","UnlA_8")
LabArep2E=GetExprssionData(LabArep2)[,c(1,2,3,7,8,9)]
colnames(LabArep2E)=c("LabA_0_2","LabA_4_2","LabA_8_2","UnlA_0_2","UnlA_4_2","UnlA_8_2")

LabCrep1E=GetExprssionData(LabCrep1)[,c(1,2,3,7,8,9)]
colnames(LabCrep1E)=c("LabC_0","LabCA_4","LabC_8","UnlC_0","UnlC_4","UnlC_8")
LabCrep2E=GetExprssionData(LabCrep2)[,c(1,2,3,7,8,9)]
colnames(LabCrep2E)=c("LabC_0_2","LabC_4_2","LabC_8_2","UnlC_0_2","UnlC_4_2","UnlC_8_2")

LabArep1fE=ConvertENSGtoGenenameNoRownames(convertExpressToFoldCHANGE(LabArep1E)[,c(2,3,5,6)])
LabArep2fE=ConvertENSGtoGenenameNoRownames(convertExpressToFoldCHANGE(LabArep2E)[,c(2,3,5,6)])
LabCrep1fE=ConvertENSGtoGenenameNoRownames(convertExpressToFoldCHANGE(LabCrep1E)[,c(2,3,5,6)])
LabCrep2fE=ConvertENSGtoGenenameNoRownames(convertExpressToFoldCHANGE(LabCrep2E)[,c(2,3,5,6)])

z=as.data.frame(subset)
for(name in c("LabArep1fE","LabArep2fE","LabCrep1fE","LabCrep2fE")){
      z=merge(z,get(name),by="id")
}
z[,2:17]=log2(z[,2:17])
write.table(z, file = "D:/ivan/ivan6-cell/figureResult/alex/limma_subset_geneexpression/inspect_subset_geneexpression.txt", append = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")