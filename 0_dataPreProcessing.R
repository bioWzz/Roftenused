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



