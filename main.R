rm(list=ls())
library(GenomicFeatures)
library(INSPEcT)
library(Rsamtools)
library(deSolve)
library(BiocParallel)
library(compiler)
library(sqldf)
load("D:/ivan/ivan6-cell/resultOfRdata/labAandC.RData")
load("D:/ivan/ivan6-cell/resultOfRdata/HarmExpression.RData")
# mydata
# Gene
# txdb3 <- makeTxDbFromGFF(file="D:/ivan/Harm/NewTestOne/NewTest/Newgencode.v23.PR8.gtf",format="gtf",taxonomyId=9606)
# load("D:/labAandC.RData")
# LabA
# LabArep1=GetResultInspect("LabA0hpi.bam",'UnlA0hpi.bam','LabA4hpi.bam','UnlA4hpi.bam','LabA8hpi.bam','UnlA8hpi.bam',"Rep1","gene")
# LabArep2=GetResultInspect("LabA0.bam",'ULabA0.bam','LabA4.bam','ULabA4.bam','LabA8.bam','ULabA8.bam',"Rep2","gene")
# LabArep12=GetResultInspectInteger("LabA0hpi.bam",'UnlA0hpi.bam','LabA4hpi.bam','UnlA4hpi.bam','LabA8hpi.bam','UnlA8hpi.bam',"LabA0.bam",'ULabA0.bam','LabA4.bam','ULabA4.bam','LabA8.bam','ULabA8.bam',"Rep1","Rep2","gene")

# LabC
# LabCrep1=GetResultInspect("LabC0hpi.bam",'UnlC0hpi.bam','LabC4hpi.bam','UnlC4hpi.bam','LabC8hpi.bam','UnlC8hpi.bam',"Rep1","gene")
# LabCrep2=GetResultInspect("LabC0.bam",'ULabC0.bam','LabC4.bam','ULabC4.bam','LabC8.bam','ULabC8.bam',"Rep2","gene")
# LabCrep12=GetResultInspectInteger("LabC0hpi.bam",'UnlC0hpi.bam','LabC4hpi.bam','UnlC4hpi.bam','LabC8hpi.bam','UnlC8hpi.bam',"LabC0.bam",'ULabC0.bam','LabC4.bam','ULabC4.bam','LabC8.bam','ULabC8.bam',"Rep1","Rep2","gene")
compareofInspectResult(LabArep1,LabArep2,"D:/ivan/ivan6-cell/figureResult/LabA")
compareofInspectResult(LabCrep1,LabCrep2,"D:/ivan/ivan6-cell/figureResult/LabC")
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#del no use data and outlier data
LabArep11=DeleteundertholdAndNAdata(ConvertFactortoNumWithRowname(LabArep1),1,9,0)
LabArep21=DeleteundertholdAndNAdata(ConvertFactortoNumWithRowname(LabArep2),1,9,0)
LabArepall=delOutliersValue(MergeByRownametwodataNoid(LabArep11,LabArep21),21,8)
LabArep12=ConvertENSGtoGenename(LabArepall[,1:21])
LabArep22=ConvertENSGtoGenename(LabArepall[,22:42])

compareofInspectResult(LabArep12,LabArep22,"D:/ivan/ivan6-cell/figureResult/delnousedata/LabA")
DawingIntegerVENNfigure(LabArep12,LabArep22,"D:/ivan/ivan6-cell/figureResult/delnousedata/LabA")
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# K-MEAN cluster
# LabArep1Kean=Kmeansclusterknum(produceKmeanData(LabArep12),0.2,10)
# LabArep2Kean=Kmeansclusterknum(produceKmeanData(LabArep22),0.2,10)
LabArep1Kean=Kmeansclusterknum(LabArep12,0.2,30,"D:/ivan/ivan6-cell/figureResult/cluster/LabA/","LabArep1Kean")
LabArep2Kean=Kmeansclusterknum(LabArep22,0.2,30,"D:/ivan/ivan6-cell/figureResult/cluster/LabA/","LabArep2Kean")

labAkmeanCompareall=allcompareBetweenKcluster(LabArep1Kean,LabArep2Kean,30,"D:/ivan/ivan6-cell/figureResult/cluster/LabA/","LabAClustercompareKean")

LabArep1K5=produceKmeanClusterResult(LabArep1Kean,6,LabArep12)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#del no use data and outlier data,将rate中na和0值删除
LabCrep11=DeleteundertholdAndNAdata(ConvertFactortoNumWithRowname(LabCrep1),1,9,0)
LabCrep21=DeleteundertholdAndNAdata(ConvertFactortoNumWithRowname(LabCrep2),1,9,0)
LabCrepall=delOutliersValue(MergeByRownametwodataNoid(LabCrep11,LabCrep21),21,8)
LabCrep12=ConvertENSGtoGenename(LabCrepall[,1:21])
LabCrep22=ConvertENSGtoGenename(LabCrepall[,22:42])

compareofInspectResult(LabCrep12,LabCrep22,"D:/ivan/ivan6-cell/figureResult/delnousedata/LabC")
DawingIntegerVENNfigure(LabCrep12,LabCrep22,"D:/ivan/ivan6-cell/figureResult/delnousedata/LabC")
# # # # # # # # # # # # # # # # # # # # 
# K-MEAN cluster
LabCrep1Kean=Kmeansclusterknum(LabCrep12,0.2,30,"D:/ivan/ivan6-cell/figureResult/cluster/LabC/","LabCrep1Kean")
LabCrep2Kean=Kmeansclusterknum(LabCrep22,0.2,30,"D:/ivan/ivan6-cell/figureResult/cluster/LabC/","LabCrep2Kean")

labCkmeanCompareall=allcompareBetweenKcluster(LabCrep1Kean,LabCrep2Kean,30,"D:/ivan/ivan6-cell/figureResult/cluster/LabC/","LabCClustercompareKean")
ResultofEachCluster(LabCrep1Kean,30,LabCrep12,path="D:/ivan/ivan6-cell/figureResult/cluster/LabC/",folder1="LabCrep1Kean")
# LabCrep1Kean4=produceKmeanClusterResult(LabCrep1Kean,6,LabCrep1)
# DrawingBoxplotandlinecartForEachClusterlog(LabCrep1Kean4,1,"D:/ivan/ivan6-cell/figureResult/cluster/LabC/","LabCrep1Kean4")


