# 出生k均值聚类的数据
produceKmeanData<-function(Expression)
{ Expression=ConvertFactortoNumWithRowname(Expression)
  Expression[,16]=Expression[,16]-Expression[,19]
  Expression[,17]=Expression[,17]-Expression[,20]
  Expression[,18]=Expression[,18]-Expression[,21]
  Expression=Expression[,c(1,2,3,4,5,6,7,8,9,16,17,18)]
}
# 对数据进行k均值聚类
# data 要聚类的数据
# rate 取数据的百分之多少进行聚类
# ratio 最小的类不少于整体的个数
# express=LabArep12
# rate=0.2
# knum=10
# path="D:/ivan/ivan6-cell/figureResult/cluster/LabA/"
# savefilename="LabArep1Kean"
Kmeansclusterknum<-function(express,rate,knum,path,savefilename) #data is a matrix contains gene features
{
  data=ConvertFactortoNumWithRowname(express)
  data=SubstantiallyGene(data,rate)
  data=ConvertFactortoNumWithRowname(data[,1:(length(data[1,])-1)])
  data=convertExpressToFoldCHANGE(data)
  n.row<-dim(data)[1]     #the number of data rows.
  result<-c()             #the results will be returned.
  for(k in 2:knum)   #number of cluster centres must lie between 1 and nrow(x),(2:nrow(x)-1)
  { 
    temp<-kmeans(data,centers=k)            #conduct k-means base on current k number.

    result<-cbind(result,temp$cluster)    #save cluster-labels
  }
  colnames(result)<-2:(dim(result)[2]+1)
  rownames(result)<-rownames(data)
  write.table(result,file=paste(path,savefilename,knum,".txt",sep=""),sep = "\t",row.names = TRUE,col.names = TRUE)
  return(result)
}
# 对数据进行k均值聚类
# data 要聚类的数据
# rate 取数据的百分之多少进行聚类
# ratio 最小的类不少于整体的个数
# express=LabArep12
# rate=0.2
# ratio=0.02 
Kmeanscluster<-function(express,rate,ratio) #data is a matrix contains gene features
{
  data=ConvertFactortoNumWithRowname(express)
  data=SubstantiallyGene(data,rate)
  data=ConvertFactortoNumWithRowname(data[,1:(length(data[1,])-1)])
  data=convertExpressToFoldCHANGE(data)
  n.row<-dim(data)[1]     #the number of data rows.
  result<-c()             #the results will be returned.
  for(k in 2:(n.row-1))   #number of cluster centres must lie between 1 and nrow(x),(2:nrow(x)-1)
  { 
    temp<-kmeans(data,centers=k)            #conduct k-means base on current k number.
    min.cluster.size<-min(temp$size)        #minimum cluster size
    if((min.cluster.size>3))       #criteria: minimum cluster size accounts for more than 2 percent of the total.
    {
      result<-cbind(result,temp$cluster)    #save cluster-labels
    }
    else{break}
  }
  colnames(result)<-2:(dim(result)[2]+1)
  rownames(result)<-rownames(data)
  return(result)
}
# 提取k均值聚类的结果
# kresult 聚类的类别结果
# clusterNum 累的个数
# data 聚类原始数据
produceKmeanClusterResult<-function(kresult,clusterNum,data)
{
  kresult=as.data.frame(kresult[,clusterNum-1])
  rownames(kresult)<-rownames(kresult)
  colnames(kresult)<-"KmeanCluster"
  data1=MergeByRownametwodataNoid(kresult,data) 
  return(data1)
}
# rep1和rep2聚类结果的比较
# KclustResult1=LabArep1Kean
# KclustResult2=LabArep2Kean
# k=2
compareBetweenKcluster<-function(KclustResult1,KclustResult2,k)
{ colnum=length(KclustResult1[1,])
  all=MergeByRownametwodataNoid(as.data.frame(KclustResult1),as.data.frame(KclustResult2))
  KclustResult1=all[,1:colnum]
  KclustResult2=all[,(colnum+1):(colnum*2)]
  
  mat=matrix(nrow=k,ncol=k) 
  dat=as.data.frame(mat) 
  for(i in 1:k){
    colnames(dat)[i]=paste("rep1cluster",i,sep="")
    rownames(dat)[i]=paste("rep2cluster",i,sep="")
  }
  for(i in 1:k){
    
    index<-KclustResult1[,k-1]==i
    s1=as.data.frame(KclustResult1[index,])
    colnames(dat)[i]=paste(colnames(dat)[i],"_",as.character(length(s1[,1])),sep="")
    for(j in 1:k){
      
      index1<-KclustResult2[,k-1]==j
      s2=as.data.frame(KclustResult2[index1,])
      rep1andrep2=intersect(as.vector(rownames(s1)),as.vector(rownames(s2)))
      if(i==1){rownames(dat)[j]=paste(rownames(dat)[j],"_",as.character(length(s2[,1])),sep="")}
      dat[j,i]=length(rep1andrep2)
    }
  }
  return(dat)
}
allcompareBetweenKcluster<-function(KclustResult1,KcustResult2,k,path,savefilename)
{
  z=list()
  for(i in 2:k)
  {
    r=compareBetweenKcluster(KclustResult1,KcustResult2,i)
    write.table(r,file=paste(path,savefilename,k,".txt",sep=""),sep = "\t",append = TRUE,row.names = TRUE,col.names = TRUE)
    z = list(z,r)
  }
  return(z)
} 
# Drawing the boxplot figrue
# 每类的综合情况
# 提取各类
# clustersingcol,聚类标签所在的列
# DrawingBoxplotExpressionandRateUnderAllCluster
# Express=LabCrep1Kean4
# clustersingcol=1
DrawingBoxplotForEachClusterlog<-function(Express,clustersingcol,savefilename)
{
  clutser=unique(unlist(Express[,clustersingcol]))
  z=ceiling(sqrt(length(clutser))) 
  graphics.off()
  # png(file=savefilename, bg = "white")
  par(mfrow=c(z,z))
  for(i in clutser){
    index<-Express[,clustersingcol]==i
    assign(paste("data",i,sep=""),Express[index,2:(length(Express[1,]))])
    boxplot(log2(get(paste("data",i,sep=""))),col=c("steelblue"),ylim = c(-10, 10),ylab=paste("Cluster",i,sep=""),xlab="synthesis  degration  processing  4su_total  4su_pre  total_total  total_pre",las=1, font.lab=0.5)
  }
  # dev.off()
 }
DrawingBoxplotForEachCluster<-function(Express,clustersingcol,savefilename)
{
  clutser=unique(unlist(Express[,clustersingcol]))
  z=ceiling(sqrt(length(clutser))) 
  graphics.off()
  # png(file=savefilename, bg = "white")
  par(mfrow=c(z,z))
  for(i in clutser){
    index<-Express[,clustersingcol]==i
    assign(paste("data",i,sep=""),Express[index,2:(length(Express[1,]))])
    boxplot((get(paste("data",i,sep=""))),col=c("steelblue"),ylim = c(-10, 10),ylab=paste("Cluster",i,sep=""),xlab="synthesis  degration  processing  4su_total  4su_pre  total_total  total_pre",las=1, font.lab=0.5)
  }
  # dev.off()
}
# foldername=folder

DrawingBoxplotandlinecartForEachClusterlog<-function(Express,clustersingcol,k,path,foldername)
{  
   clutser=unique(unlist(Express[,clustersingcol]))
   for(i in clutser){
    graphics.off()
    plot.new()
    par(mar=c(1,1,1,1))
    png(file=paste(path,foldername,"/",i,".png",sep=""), bg = "white",width = 900, height =480)
    index<-Express[,clustersingcol]==i
    assign(paste("data",i,sep=""),Express[index,2:(length(Express[1,]))])
    tep=get(paste("data",i,sep=""))
    par(fig=c(0,1,0,1))
    boxplot(log2(get(paste("data",i,sep=""))),col=c("steelblue"),ylim = c(-10, 10),ylab=paste("Cluster",i,sep=""),xlab="synthesis  degration  processing  4su_total  4su_pre  total_total  total_pre",xaxt='n',font.lab=3)
    if(length(tep[,1])>3){
    for(j in seq(1,length(tep[1,]),3))
    { 
      x=colMeans(as.matrix(sapply(tep[,j:(j+2)],as.numeric)))
      y=c(1,4,8)
      par(fig=c(0.01+(0.14*(as.integer(j/3))),0.15+(0.14*(as.integer(j/3))),0,0.3), new=TRUE)
      plot(y, x, type = "n",xaxt='n',yaxt='n',ann=FALSE)
      points(y, x, type = "o", pch = 1, col = "red", lty = 1, lwd = 2)
    }
    }
   dev.off()
  }
}

# 生成所有的聚类结果和图


# Keancluster=LabCrep1Kean
# kmun=30
# data=LabCrep12
# path="D:/ivan/ivan6-cell/figureResult/cluster/LabC/"
# folder1="LabCrep1Kean"

ResultofEachCluster<-function(Keancluster,kmun,data,path,folder1){
  for(k in 2:kmun){
    folder=paste(folder1,k,sep="")
    dir.create(paste(path,folder,"/",sep=""))
    singleclusterR=produceKmeanClusterResult(Keancluster,k,data)
    write.table(singleclusterR,file=paste(path,folder,"/Cluster",k,".txt",sep=""),sep = "\t",append =F,row.names = TRUE,col.names = TRUE)
    DrawingBoxplotandlinecartForEachClusterlog(singleclusterR,1,k,path,folder)
  }
}



