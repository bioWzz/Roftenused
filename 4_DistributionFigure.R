# 在一个画布上画多张图，数据进行了log2转化
ManyDistributionFigurelog(Express,"red","Fraction of Gene",0,3,0,3)
Express=as.data.frame(LabArep1[,1:9])
savefilename="D:/ivan/ivan6-cell/figureResult/LabA"
linecolour="red"
mynrow=3
myncol=3
myxlab="Fraction of Gene"
minxlim=0
maxxlim=3
minylim=0
maxylim=3
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


