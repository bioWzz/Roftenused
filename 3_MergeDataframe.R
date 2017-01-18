# 删除多余变量
# del="^p"
delvari<-function(del){
vari=objects()
delva=as.vector(grep(del,vari,invert=TRUE))  
vari=vari[delva]
rm(list=vari)
}

#Get Expression Data
GetExprssionData<-function(Express){
  Express=Express[,10:length(Express[1,])]
  return(Express)
}

# all the weight to raw data
ChangeWeightsOfExpression<-function(Express,weigth)
{ 
  Express=ConvertFactortoNumWithRowname(Express)
  for(i in 1:length(Express[1,])){
    Express[,i]=Express[,i]*weigth[i]
  }
  return(Express)
}

#produce the new inpsect data
produceNewRate<-function(Express)
{
  Express=ConvertFactortoNumWithRowname(Express)
  foursu_exons=Express[,1:3]
  total_exons=Express[,7:9]
  foursu_introns=Express[,4:6]
  total_introns=Express[,10:12]
  tpts <- c(0,4,8)
  tL <- 4
  
  result<-newINSPEcT(tpts, tL, foursu_exons, total_exons,
                     foursu_introns, total_introns, BPPARAM=SerialParam())
  Express=result@ratesFirstGuess@assayData$exprs
  Express=Express[,7:length(Express[1,])]
  
  Express1=MergeByRownametwodata(Express,foursu_exons,"id")
  Express1=MergeByRownametwodata(Express1,foursu_introns,"id")
  Express1=MergeByRownametwodata(Express1,total_exons,"id")
  Express1=MergeByRownametwodata(Express1,total_introns,"id")
  return(Express1)
}

# 为rep1，rep2画散点图及分布图
compareofInspectResult<-function(Express1,Express2,filename)
{
  # # compare Expression
  DrawingIntegerScatterPlot(Express1,Express2,21,10,21,paste(filename,"Expression.png",sep=""),4,3,"ix",c(0,1000),c(0,1000),"green")
  ManyDistributionFigurelog(as.data.frame(Express1[,10:21]),paste(filename,"Rep1ExpressionDistribution.png",sep=""),"red",4,3,"Fraction of Gene",0,15,0,1)
  ManyDistributionFigurelog(as.data.frame(Express2[,10:21]),paste(filename,"Rep2ExpressionDistribution.png",sep=""),"red",4,3,"Fraction of Gene",0,15,0,1)
  # compare Rate
  DrawingIntegerScatterPlot(Express1,Express2,21,1,9,paste(filename,"Rate.png",sep=""),3,3,"ix",c(0,200),c(0,200),"green")
  ManyDistributionFigurelog(as.data.frame(Express1[,1:9]),paste(filename,"Rep1RateDistribution.png",sep=""),"red",3,3,"Fraction of Gene",0,10,0,1)
  ManyDistributionFigurelog(as.data.frame(Express2[,1:9]),paste(filename,"Rep2RateDistribution.png",sep=""),"red",3,3,"Fraction of Gene",0,10,0,1)
  
}

# get the values which are same to harm
getvaluelog2<-function(Express){
  Express=log2(ConvertFactortoNumWithRowname(Express[,c(10,11,12,16,17,18)]))
  return(Express)
}
getvaluelog2subset<-function(Express,subset){
  Express=log2(ConvertFactortoNumWithRowname(Express[,c(10,11,12,16,17,18)]))
  Result=ConvertENSGtoGenenameNoRownames(Express)
  colnames(subset)="id"
  subsetResult=merge(subset,Result,by="id")
  rownames(subsetResult)=subsetResult[,1]
  subsetResult=subsetResult[,2:(length(subsetResult[1,]))]
  return(subsetResult)
}

# hebing wan fen kai fan hui
mergeandchooseOneData<-function(Express1,Express2,which){
  colnum=length(Express1[1,])
  all=MergeByRownametwodata(as.data.frame(Express1),as.data.frame(Express2),"id")
  one=all[,1:colnum]
  two=all[,(colnum+1):(colnum*2)]
  if(which==1){
    return(one)
  }
  if(which==2){
    return(two)
  }
}

# Merge of many dataframe
# Merge dataframes which have the same lines,Named row name as the first colcume
# 把数据框结合起来，并将行名提取出来作为第一列，并给第一列起id
# NameofMerge=c("data1","data2",...) 数据框的个数大于等于2个
# FirstColnamesId 第一列的列名
simpleMergeWFistCol<-function(NameofMerge,FirstColnamesId)
{
  result=get(NameofMerge[1])
  
  for(i in 2:length(NameofMerge)){
     result=cbind(result,get(NameofMerge[i]))
  }
  result=cbind(rownames(result),result)
  colnames(result)[1]=FirstColnamesId
  return(result)
}
# Merge dataframes which have the same lines
# 把数据框结合起来
# NameofMerge=c("data1","data2",...) 数据框的个数大于等于2个
simpleMerge<-function(NameofMerge)
{
  result=get(NameofMerge[1])
  
  for(i in 2:length(NameofMerge)){
    result=cbind(result,get(NameofMerge[i]))
  }
  return(result)
}

# Merge by rownames
# 把数据框合并起来按照行名，返回合并起来的数据，并且行名作为第一列
MergeByRowname<-function(NameofMerge,FirstColnamesId)
{ 
  result=get(NameofMerge[1])
  result=cbind(rownames(result),result)
  colnames(result)[1]=FirstColnamesId
  for(i in 2:length(NameofMerge)){
    z=get(NameofMerge[i])
    z=cbind(rownames(z),z)
    colnames(z)[1]=FirstColnamesId
   result=merge(result,z,by=FirstColnamesId)
  }
}

# 把数据框合并起来按照行名，返回合并起来的数据,数据有行名
# NameofMerge=c("Express","foursu_exons","foursu_introns","total_exons","total_introns")
# FirstColnamesId="id"
MergeByRowname<-function(NameofMerge,FirstColnamesId)
{ 
  result=get(NameofMerge[1])
  result=cbind(rownames(result),result)
  colnames(result)[1]=FirstColnamesId
  for(i in 2:length(NameofMerge)){
    z=get(NameofMerge[i])
    z=cbind(rownames(z),z)
    colnames(z)[1]=FirstColnamesId
    result=merge(result,z,by=FirstColnamesId)
  }
  rownames(result)=result[,1]
  result=result[,2:length(result[1,])]
  
  return(result)
}

# 两个数据集合合并
MergeByRownametwodata<-function(data1,data2,FirstColnamesId)
{ 

  data1=cbind(rownames(data1),data1)
  colnames(data1)[1]=FirstColnamesId
 
  data2=cbind(rownames(data2),data2)
  colnames(data2)[1]=FirstColnamesId

  result=merge(data1,data2,by=FirstColnamesId)
  rownames(result)=result[,1]
  result=result[,2:length(result[1,])]
  return(result)
}

# 两个数据集合合并,不需要id
MergeByRownametwodataNoid<-function(data1,data2)
{ 
  
  data1=cbind(rownames(data1),data1)
  colnames(data1)[1]="id"
  
  data2=cbind(rownames(data2),data2)
  colnames(data2)[1]="id"
  
  result=merge(data1,data2,by="id")
  rownames(result)=result[,1]
  result=result[,2:length(result[1,])]
  return(result)
}

# 删除小于阈值的数据,
DeleteundertholdAndNAdata<-function(Express,Scol,Ecol,value){
  Express=na.omit(Express)
  Express[apply(Express[,Scol:Ecol]<=value, FUN = any, 1), ] = NA
  Express=na.omit(Express)
  return(Express)
}

# 删除小于0的数据,
DeleteunderZeroAndNAdata<-function(Express){
  Express=na.omit(Express)
  Express[apply(Express[,1:length(Express[1,])]<=0, FUN = any, 1), ] = NA
  Express=na.omit(Express)
  return(Express)
}

# convert the factor data to num data 
ConvertFactortoNum<-function(Express)
{ 
  if(length(Express[1,])==1){
    result=as.data.frame(as.numeric(as.character(Express[,1])))
    result=as.data.frame(result)
  }else{
    result=as.data.frame(as.numeric(as.character(Express[,1])))
    for(i in 2:length(Express[1,]))
     {
       result=cbind(result,as.data.frame(as.numeric(as.character(Express[,i]))))
     }
    result=as.data.frame(result)
  }
  return(result)
}

# convert the factor data to num data with rowname 
ConvertFactortoNumWithRowname<-function(Express)
{ 
  if(length(Express[1,])==1){
    result=as.data.frame(as.numeric(as.character(Express[,1])))
    result=as.data.frame(result)
    rownames(result)=rownames(Express)
    colnames(result)=colnames(Express)
  }else{
  result=as.data.frame(as.numeric(as.character(Express[,1])))
  for(i in 2:length(Express[1,]))
  {
    result=cbind(result,as.data.frame(as.numeric(as.character(Express[,i]))))
  }
  result=as.data.frame(result)
  rownames(result)=rownames(Express)
  colnames(result)=colnames(Express)}
  result=as.data.frame(result)
  return(result)
}


# convert ENSG TO GeneName
# p=ConvertENSGtoGenename(LabArep12)
ConvertENSGtoGenename<-function(Express)
{
  ENSG_GeneNames <- read.delim("D:/ivan/ivan6-cell/figureResult/ENSG_GeneNames.txt")
  colnames(ENSG_GeneNames)[1]="id"
  name=as.data.frame(rownames(Express))
  name1=as.data.frame(strsplit(as.character(name[,1]),'\\.'))
  name1<- t(name1)
  name1=as.data.frame(name1[,1])
  colnames(name1)[1]="id"
  Express=cbind(name1,Express)
  ALL=merge(Express,ENSG_GeneNames,by="id")
  rownames(ALL)=ALL[,25]
  ALL=ALL[,2:22]
}

#gene symbol zuowei diyilie 
ConvertENSGtoGenenameNoRownames<-function(Express)
{
  ENSG_GeneNames <- read.delim("D:/ivan/ivan6-cell/figureResult/ENSG_GeneNames.txt")
  colnames(ENSG_GeneNames)[1]="id"
  name=as.data.frame(rownames(Express))
  name1=as.data.frame(strsplit(as.character(name[,1]),'\\.'))
  name1<- t(name1)
  name1=as.data.frame(name1[,1])
  colnames(name1)[1]="id"
  Express=cbind(name1,Express)
  ALL=merge(Express,ENSG_GeneNames,by="id")
  ALL[,1]=ALL[,length(ALL[1,])]
  ALL=ALL[,1:(length(ALL[1,])-3)]
}

# delete the no use data of dataframe
# threshold value dataframe 中值应该满足的条件
# deltype： all 所有的 ，any任何一个
delNoUseData<-function(Express,Scol,Ecol,threshold,deltype)
{
  Express=na.omit(Express)
  Express=ConvertFactortoNumWithRowname(Express)
  Express[apply(Express[,Scol:Ecol]<=threshold, FUN = deltype, 1), ] = NA
  Express=na.omit(Express)
}

# 删除异常值
delOutliersValue<-function(data,Numofcols,Threshold){
  library(matrixStats)
  library(Biobase)
  data=ConvertFactortoNumWithRowname(data)
  data=na.omit(data)
  for(i in 1:Numofcols){
  data=na.omit(data)
  MyMatrix=as.matrix(cbind(data[,i],data[,(Numofcols+i)]))
  min=as.data.frame(rowMin(MyMatrix))
  max=as.data.frame(rowMax(MyMatrix))
  index<-max==0
  max[index]=0.0001
  min[index]=0.0001
  index.enst<-((min/max)<(1/Threshold))    #the index of enst for one ensg
  data[index.enst,]=NA
  }
  data=na.omit(data)
  return(data)
}

# 提取高表达基因
SubstantiallyGene<-function(Express,rate)
{
  Express=ConvertFactortoNumWithRowname(Express)
  Express=na.omit(Express)
  Express[apply(Express<=0, FUN = any, 1), ] = NA
  Express=na.omit(Express)
  allExpression=as.data.frame(rowSums(Express[,10:(length(Express[1,]))]))
  Express=cbind(as.data.frame(Express),as.data.frame(allExpression))
  Express=Express[order(Express[,length(Express[1,])],decreasing=TRUE),]
  highexpressiongene=Express[1:as.integer(length(Express[,1])*rate),]
  return(highexpressiongene)
}

# 画韦恩图
# Express1=LabArep1
# Express2=LabArep2
DawingVENNfiguroftwoset<-function(Express1,Express2){
  library(gplots)
  i=1
  GroupA<-as.list(rownames(Express1)) 
  GroupB<-as.list(rownames(Express2))
  input <-list(GroupA,GroupB)
  venn(input)
  
}

#两组数据集递进画韦恩图
# # 画韦恩图
# # 画韦恩图,得到rep1和rep2在按照表达从高到低排序的，取相同比例的子集中的交集的韦恩图
DawingIntegerVENNfigure<-function(Express1,Express2,savefilename){
graphics.off()
png(file=paste(savefilename,"VENN.png",sep=""), bg = "white")
par(mfrow=c(3,3))
for(i in 1:9){
  highLabArep12=SubstantiallyGene(Express1,0.1*i)
  highLabArep22=SubstantiallyGene(Express2,0.1*i)
  DawingVENNfiguroftwoset(highLabArep12,highLabArep22)
 }
dev.off()
}

# 将表达水平编程fold倍数，和0时刻相比
# Express=subLabA_REP1
convertExpressToFoldCHANGE<-function(Express){
  Express=ConvertFactortoNumWithRowname(Express)
  Express[Express==0]<-0.0001
  type=c()
  for(i in seq(1,length(Express[1,]),3)){
    type=c(type,i)
  }
  for(i in type){
    Express[,i+1]=log2(Express[,i+1]/Express[,i])
    Express[,i+2]=log2(Express[,i+2]/Express[,i])
    Express[,i]=0
  }
  return(Express)
}

convertExpressToFoldCHANGElog<-function(Express){
  Express=ConvertFactortoNumWithRowname(Express)
  type=c()
  for(i in seq(1,length(Express[1,]),3)){
    type=c(type,i)
  }
  for(i in type){
    Express[,i+1]=2^(Express[,i+1] - Express[,i])
    Express[,i+2]=2^(Express[,i+2] - Express[,i])
    Express[,i]=0
  }
  return(Express)
}

ConverttoDataWithRownames<-function(Express){
  rownames(Express)=Express[,1]
  Express=Express[,2:length(Express[1,])]
  return(Express)
}

