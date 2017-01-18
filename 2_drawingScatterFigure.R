# drawing the scatter plot between two data frames of data1 and data2
# data1有行名，列是数据
# data2有行名，列是数据
# rep1第一个数据,一列数据，有行名
# rep2第二个数据，一列数据，有行名
# lines 画布的行数
# colnames画布的列数
# title x轴下方的标记
# xlim 坐标范围
# graphics.off()
# par(mfrow=c(2,3))
# rep1=as.data.frame(rep1foursu_exons)
# rep2=as.data.frame(rep1foursu_exons)
# DrawingSingleScatterPlot(rep1,rep2,"IT",c(0,50),c(0,50),"red")

# lengthofdata1数据框的长度
#SRawingnumofcols 从那列开始画
#ERawingnumofcols 画到那列结束
# 两个数据框画多个散点图
# data1=LabArep1E
# data2=LabArep2E
# lengthofdata1=12
# SRawingnumofcols=1
# ERawingnumofcols=12
# rowlines=4
# collines=3
# title="ix"
# myxlim=c(0,1000)
# myylim=c(0,1000)
# linecolour="green"
# savefilename="LabA"

DrawingIntegerScatterPlot<-function(data1,data2,lengthofdata1,SRawingnumofcols,ERawingnumofcols,savefilename,rowlines,collines,title,myxlim,myylim,linecolour)
{ 
  graphics.off()
  png(file=savefilename, bg = "white")
  par(mfrow=c(rowlines,collines))
  data1=ConvertFactortoNumWithRowname(data1)
  data2=ConvertFactortoNumWithRowname(data2)
  testl=MergeByRownametwodata(data1,data2,"id")
  testlabArep1=testl[,1:lengthofdata1]
  testlabArep2=testl[,(lengthofdata1+1):(lengthofdata1*2)]
  data1=as.data.frame(testlabArep1[,SRawingnumofcols:ERawingnumofcols])
  data2=as.data.frame(testlabArep2[,SRawingnumofcols:ERawingnumofcols])
  
  for(i in 1:length(data1[1,]))
  { print(i)
    rep1=ConvertFactortoNumWithRowname(as.data.frame(data1[,i]))
    rownames(rep1)=rownames(data1)
    rep2=ConvertFactortoNumWithRowname(as.data.frame(data2[,i]))
    rownames(rep2)=rownames(data2)
    DrawingSingleScatterPlot(rep1,rep2,title,myxlim,myylim,linecolour,as.character(colnames(data1)[i]),as.character(colnames(data2)[i]))
  }
  dev.off()
}
# 两列数画一个散点图
DrawingSingleScatterPlot<-function(rep1,rep2,title,myxlim,myylim,linecolour,myxlab,myylab)
{
  library(ggplot2)
  rep1=ConvertFactortoNumWithRowname(rep1)
  rep2=ConvertFactortoNumWithRowname(rep2)
  rep1=cbind(rownames(rep1),rep1)
  colnames(rep1)[1]="id"
  
  rep2=cbind(rownames(rep2),rep2)
  colnames(rep2)[1]="id"
  
  data=as.data.frame(merge(rep1,rep2, by ="id"))
  data=data[,c(2,3)]
  data=ConvertFactortoNum(data)
  
  data=DeleteunderZeroAndNAdata(data)
  model<-lm(data[,2]~data[,1])
  summary(model)
  model$coefficients[[2]]
  plot(data[,1],data[,2],xlim=myxlim,ylim=myylim,main=as.character( model$coefficients[[2]]),ylab=myylab, xlab=myxlab)
  abline(model,col=linecolour)
}

# 在一个画布上画多张图，数据进行了log2转化
# ManyDistributionFigurelog(Express,"red","Fraction of Gene",0,3,0,3)
# Express=as.data.frame(LabArep1[,1:9])
# savefilename="D:/ivan/ivan6-cell/figureResult/LabA"
# linecolour="red"
# mynrow=3
# myncol=3
# myxlab="Fraction of Gene"
# minxlim=0
# maxxlim=3
# minylim=0
# maxylim=3
ManyDistributionFigurelog<-function(Express,savefilename,linecolour,mynrow,myncol,myxlab,minxlim,maxxlim,minylim,maxylim){
  Express=ConvertFactortoNumWithRowname(Express)
  colname=colnames(Express)
  library(ggplot2)
  library(cowplot)
  graphics.off()
  
  plots <- NULL
  for(i in 1:length(Express[1,]))
  { 
    tmp=as.data.frame(Express[,i])
    e1=ConvertFactortoNum(tmp)
    e1=as.data.frame(as.numeric(unlist(e1))) 
    e1=e1[e1>0]
    e1=log2(e1)
    data <- data.frame(x = e1)
    plots[[i]]<- ggplot(data=data, aes(x = x, y = ..density..))+geom_density(data=data,colour = linecolour,size = 0.3)+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x = element_text(size =5,face = "bold"),axis.title.y = element_text(size = 5,face = "bold"),axis.text.x = element_text(size = 5,face = "bold"),axis.text.y = element_text(size = 5,face = "bold"))+ xlab(as.character(colnames(Express)[i])) + ylab(myxlab)+xlim(minxlim,maxxlim)+ylim(minylim,maxylim)
  }
  
  ggsave(plot_grid(plotlist = plots,nrow =mynrow,ncol = myncol), file=savefilename, width=4, height=4)
}

# 在一个画布上画多张图，数据不进行log2转化
ManyDistributionFigure<-function(Express,linecolour,mynrow,myncol,myxlab,minxlim,maxxlim,minylim,maxylim){
  Express=ConvertFactortoNumWithRowname(Express)
  colname=colnames(Express)
  library(ggplot2)
  library(cowplot)
  graphics.off()
  plots <- NULL
  for(i in 1:length(Express[1,]))
  { 
    tmp=as.data.frame(Express[,i])
    e1=ConvertFactortoNum(tmp)
    e1=as.data.frame(as.numeric(unlist(e1))) 
    e1=e1[e1>0]
    data <- data.frame(x = e1)
    plots[[i]]<- ggplot(data=data, aes(x = x, y = ..density..))+geom_density(data=data,colour = linecolour,size = 1)+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))+ xlab(as.character(colnames(Express)[i])) + ylab(myxlab)+xlim(minxlim,maxxlim)+ylim(minylim,maxylim)
  }
  plot_grid(plotlist = plots,nrow =mynrow,ncol = myncol)
}


