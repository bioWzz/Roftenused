# rep1 or rep2 result

GetResultInspect<-function(suh0,total0,suh4,total4,suh8,total8,rep,genetype)
{

  # read the data
  # Lab-0
  paths_4su_0 <- system.file(paste('extdata/',rep,sep=""),suh0, package="INSPEcT")
  paths_total_0 <- system.file(paste('extdata/',rep,sep=""),total0, package="INSPEcT")
  H0hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su_0, paths_total_0,by =genetype)
  # Lab-4
  paths_4su_4 <- system.file(paste('extdata/',rep,sep=""),suh4, package="INSPEcT")
  paths_total_4 <- system.file(paste('extdata/',rep,sep=""),total4, package="INSPEcT")
  H4hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su_4, paths_total_4,by =genetype)
  # Lab-8
  paths_4su_8 <- system.file(paste('extdata/',rep,sep=""),suh8, package="INSPEcT")
  paths_total_8 <- system.file(paste('extdata/',rep,sep=""),total8, package="INSPEcT")
  H8hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su_8, paths_total_8,by =genetype)
  

  ######################
  ##LabA rate
  
  H0hpifoursu_exons=H0hpiRPKMsOut$rpkms$foursu_exons
  H4hpifoursu_exons=H4hpiRPKMsOut$rpkms$foursu_exons
  H8hpifoursu_exons=H8hpiRPKMsOut$rpkms$foursu_exons
  foursu_exons=cbind(H0hpifoursu_exons,H4hpifoursu_exons,H8hpifoursu_exons)
  colnames(foursu_exons) <-c("foursu_exonsOnehour","foursu_exonsFourhour","foursu_exonsEighthour")
  
  H0hpifoursu_introns=H0hpiRPKMsOut$rpkms$foursu_introns
  H4hpifoursu_introns=H4hpiRPKMsOut$rpkms$foursu_introns
  H8hpifoursu_introns=H8hpiRPKMsOut$rpkms$foursu_introns
  foursu_introns=cbind(H0hpifoursu_introns,H4hpifoursu_introns,H8hpifoursu_introns)
  colnames(foursu_introns) <-c("foursu_intronsOnehour","foursu_intronsFourhour","foursu_intronsEighthour")
  
  # UnlA0hpitotal_exons=Rep1LabA0hpiRPKMsOut$rpkms$total_exons
  TH0hpitotal_exons=H0hpiRPKMsOut$rpkms$total_exons
  TH4hpitotal_exons=H4hpiRPKMsOut$rpkms$total_exons
  TH8hpitotal_exons=H8hpiRPKMsOut$rpkms$total_exons
  total_exons=cbind(TH0hpitotal_exons,TH4hpitotal_exons,TH8hpitotal_exons)
  colnames(total_exons) <-c("total_exonsOnehour","total_exonsFourhour","total_exonsEighthour")
  
  TH0hpitotal_introns=H0hpiRPKMsOut$rpkms$total_introns
  TH4hpitotal_introns=H4hpiRPKMsOut$rpkms$total_introns
  TH8hpitotal_introns=H8hpiRPKMsOut$rpkms$total_introns
  total_introns=cbind(TH0hpitotal_introns,TH4hpitotal_introns,TH8hpitotal_introns)
  colnames(total_introns) <-c("total_intronsOnehour","total_intronsFourhour","total_intronsEighthour")
  
  tpts <- c(0,4,8)
  tL <- 4

  result<-newINSPEcT(tpts, tL, foursu_exons, total_exons,
                                               foursu_introns, total_introns, BPPARAM=SerialParam())
  print(paths_4su_8)
  Express=result@ratesFirstGuess@assayData$exprs
  Express=Express[,7:length(Express[1,])]
  
  Express1=MergeByRownametwodata(Express,foursu_exons,"id")
  Express1=MergeByRownametwodata(Express1,foursu_introns,"id")
  Express1=MergeByRownametwodata(Express1,total_exons,"id")
  Express1=MergeByRownametwodata(Express1,total_introns,"id")
  return(Express1)
}


GetcountResultInspect<-function(suh0,total0,suh4,total4,suh8,total8,rep,genetype)
{
  library(plyr)
  # read the data
  # Lab-0
  paths_4su_0 <- system.file(paste('extdata/',rep,sep=""),suh0, package="INSPEcT")
  paths_total_0 <- system.file(paste('extdata/',rep,sep=""),total0, package="INSPEcT")
  H0hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su_0, paths_total_0,by =genetype)
  # Lab-4
  paths_4su_4 <- system.file(paste('extdata/',rep,sep=""),suh4, package="INSPEcT")
  paths_total_4 <- system.file(paste('extdata/',rep,sep=""),total4, package="INSPEcT")
  H4hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su_4, paths_total_4,by =genetype)
  # Lab-8
  paths_4su_8 <- system.file(paste('extdata/',rep,sep=""),suh8, package="INSPEcT")
  paths_total_8 <- system.file(paste('extdata/',rep,sep=""),total8, package="INSPEcT")
  H8hpiRPKMsOut <- makeRPKMs(txdb3, paths_4su_8, paths_total_8,by =genetype)
  
  
  ######################
  ##LabA rate
  
  H0hpifoursu_exons=H0hpiRPKMsOut$counts$foursu$exonCounts
  H4hpifoursu_exons=H4hpiRPKMsOut$counts$foursu$exonCounts
  H8hpifoursu_exons=H8hpiRPKMsOut$counts$foursu$exonCounts
  foursu_exons=cbind(rownames(H0hpifoursu_exons),H0hpifoursu_exons,H4hpifoursu_exons,H8hpifoursu_exons)
  colnames(foursu_exons) <-c("id","foursu_exonsOnehour","foursu_exonsFourhour","foursu_exonsEighthour")
  
  H0hpifoursu_introns=H0hpiRPKMsOut$counts$foursu$intronCounts
  H4hpifoursu_introns=H4hpiRPKMsOut$counts$foursu$intronCounts
  H8hpifoursu_introns=H8hpiRPKMsOut$counts$foursu$intronCounts
  foursu_introns=cbind(rownames(H0hpifoursu_introns),H0hpifoursu_introns,H4hpifoursu_introns,H8hpifoursu_introns)
  colnames(foursu_introns) <-c("id","foursu_intronsOnehour","foursu_intronsFourhour","foursu_intronsEighthour")
  
  # UnlA0hpitotal_exons=Rep1LabA0hpiRPKMsOut$rpkms$total_exons
  TH0hpitotal_exons=H0hpiRPKMsOut$counts$total$exonCounts
  TH4hpitotal_exons=H4hpiRPKMsOut$counts$total$exonCounts
  TH8hpitotal_exons=H8hpiRPKMsOut$counts$total$exonCounts
  total_exons=cbind(rownames(TH0hpitotal_exons),TH0hpitotal_exons,TH4hpitotal_exons,TH8hpitotal_exons)
  colnames(total_exons) <-c("id","total_exonsOnehour","total_exonsFourhour","total_exonsEighthour")
  
  TH0hpitotal_introns=H0hpiRPKMsOut$counts$total$intronCounts
  TH4hpitotal_introns=H4hpiRPKMsOut$counts$total$intronCounts
  TH8hpitotal_introns=H8hpiRPKMsOut$counts$total$intronCounts
  total_introns=cbind(rownames(TH0hpitotal_introns),TH0hpitotal_introns,TH4hpitotal_introns,TH8hpitotal_introns)
  colnames(total_introns) <-c("id","total_intronsOnehour","total_intronsFourhour","total_intronsEighthour")
  
  restult=merge(foursu_exons, foursu_introns, by ="id",all = TRUE,incomparables = NULL)
  restult=merge(restult, total_exons, by ="id",all = TRUE,incomparables = NULL) 
  restult=merge(restult, total_introns, by ="id",all = TRUE,incomparables = NULL) 
  return(restult)
}

# rep1 and rep2 result

GetResultInspectInteger<-function(suh0,total0,suh4,total4,suh8,total8,suh01,total01,suh41,total41,suh81,total81,rep,rep1,genetype)
{
  # rep1
  # Lab-0
  rep1paths_4su_0 <- system.file(paste('extdata/',rep,sep=""),suh0, package="INSPEcT")
  rep1paths_total_0 <- system.file(paste('extdata/',rep,sep=""),total0, package="INSPEcT")
  rep1H0hpiRPKMsOut <- makeRPKMs(txdb3, rep1paths_4su_0, rep1paths_total_0,by =genetype)
  # Lab-4
  rep1paths_4su_4 <- system.file(paste('extdata/',rep,sep=""),suh4, package="INSPEcT")
  rep1paths_total_4 <- system.file(paste('extdata/',rep,sep=""),total4, package="INSPEcT")
  rep1H4hpiRPKMsOut <- makeRPKMs(txdb3, rep1paths_4su_4, rep1paths_total_4,by =genetype)
  # Lab-8
  rep1paths_4su_8 <- system.file(paste('extdata/',rep,sep=""),suh8, package="INSPEcT")
  rep1paths_total_8 <- system.file(paste('extdata/',rep,sep=""),total8, package="INSPEcT")
  rep1H8hpiRPKMsOut <- makeRPKMs(txdb3, rep1paths_4su_8, rep1paths_total_8,by =genetype)
  
  
  ######################
  ##Lab rate
  
  rep1H0hpifoursu_exons=rep1H0hpiRPKMsOut$rpkms$foursu_exons
  rep1H4hpifoursu_exons=rep1H4hpiRPKMsOut$rpkms$foursu_exons
  rep1H8hpifoursu_exons=rep1H8hpiRPKMsOut$rpkms$foursu_exons
  rep1foursu_exons=cbind(rep1H0hpifoursu_exons,rep1H4hpifoursu_exons,rep1H8hpifoursu_exons)
  colnames(rep1foursu_exons) <-c("rep1foursu_exonsOnehour","rep1foursu_exonsFourhour","rep1foursu_exonsEighthour")
  
  rep1H0hpifoursu_introns=rep1H0hpiRPKMsOut$rpkms$foursu_introns
  rep1H4hpifoursu_introns=rep1H4hpiRPKMsOut$rpkms$foursu_introns
  rep1H8hpifoursu_introns=rep1H8hpiRPKMsOut$rpkms$foursu_introns
  rep1foursu_introns=cbind(rep1H0hpifoursu_introns,rep1H4hpifoursu_introns,rep1H8hpifoursu_introns)
  colnames(rep1foursu_introns) <-c("rep1foursu_intronsOnehour","rep1foursu_intronsFourhour","rep1foursu_intronsEighthour")
  
  # UnlA0hpitotal_exons=Rep1LabA0hpiRPKMsOut$rpkms$total_exons
  rep1TH0hpitotal_exons=rep1H0hpiRPKMsOut$rpkms$total_exons
  rep1TH4hpitotal_exons=rep1H4hpiRPKMsOut$rpkms$total_exons
  rep1TH8hpitotal_exons=rep1H8hpiRPKMsOut$rpkms$total_exons
  rep1total_exons=cbind(rep1TH0hpitotal_exons,rep1TH4hpitotal_exons,rep1TH8hpitotal_exons)
  colnames(rep1total_exons) <-c("rep1total_exonsOnehour","rep1total_exonsFourhour","rep1total_exonsEighthour")
  
  rep1TH0hpitotal_introns=rep1H0hpiRPKMsOut$rpkms$total_introns
  rep1TH4hpitotal_introns=rep1H4hpiRPKMsOut$rpkms$total_introns
  rep1TH8hpitotal_introns=rep1H8hpiRPKMsOut$rpkms$total_introns
  rep1total_introns=cbind(rep1TH0hpitotal_introns,rep1TH4hpitotal_introns,rep1TH8hpitotal_introns)
  colnames(rep1total_introns) <-c("rep1total_intronsOnehour","rep1total_intronsFourhour","rep1total_intronsEighthour")
  
  # rep2
  # Lab-0
  rep2paths_4su_0 <- system.file(paste('extdata/',rep1,sep=""),suh01, package="INSPEcT")
  rep2paths_total_0 <- system.file(paste('extdata/',rep1,sep=""),total01, package="INSPEcT")
  rep2H0hpiRPKMsOut <- makeRPKMs(txdb3, rep2paths_4su_0, rep2paths_total_0,by =genetype)
  # Lab-4
  rep2paths_4su_4 <- system.file(paste('extdata/',rep1,sep=""),suh41, package="INSPEcT")
  rep2paths_total_4 <- system.file(paste('extdata/',rep1,sep=""),total41, package="INSPEcT")
  rep2H4hpiRPKMsOut <- makeRPKMs(txdb3, rep2paths_4su_4, rep2paths_total_4,by =genetype)
  # Lab-8
  rep2paths_4su_8 <- system.file(paste('extdata/',rep1,sep=""),suh81, package="INSPEcT")
  rep2paths_total_8 <- system.file(paste('extdata/',rep1,sep=""),total81, package="INSPEcT")
  rep2H8hpiRPKMsOut <- makeRPKMs(txdb3, rep2paths_4su_8, rep2paths_total_8,by =genetype)
  
  
  ######################
  ##Lab rate
  
  rep2H0hpifoursu_exons=rep2H0hpiRPKMsOut$rpkms$foursu_exons
  rep2H4hpifoursu_exons=rep2H4hpiRPKMsOut$rpkms$foursu_exons
  rep2H8hpifoursu_exons=rep2H8hpiRPKMsOut$rpkms$foursu_exons
  rep2foursu_exons1=cbind(rep2H0hpifoursu_exons,rep2H4hpifoursu_exons,rep2H8hpifoursu_exons)
  colnames(rep2foursu_exons1) <-c("rep2foursu_exonsOnehour","rep2foursu_exonsFourhour","rep2foursu_exonsEighthour")
  
  rep2H0hpifoursu_introns=rep2H0hpiRPKMsOut$rpkms$foursu_introns
  rep2H4hpifoursu_introns=rep2H4hpiRPKMsOut$rpkms$foursu_introns
  rep2H8hpifoursu_introns=rep2H8hpiRPKMsOut$rpkms$foursu_introns
  rep2foursu_introns1=cbind(rep2H0hpifoursu_introns,rep2H4hpifoursu_introns,rep2H8hpifoursu_introns)
  colnames(rep2foursu_introns1) <-c("rep2foursu_intronsOnehour","rep2foursu_intronsFourhour","rep2foursu_intronsEighthour")
  
  # UnlA0hpitotal_exons=Rep1LabA0hpiRPKMsOut$rpkms$total_exons
  rep2TH0hpitotal_exons=rep2H0hpiRPKMsOut$rpkms$total_exons
  rep2TH4hpitotal_exons=rep2H4hpiRPKMsOut$rpkms$total_exons
  rep2TH8hpitotal_exons=rep2H8hpiRPKMsOut$rpkms$total_exons
  rep2total_exons1=cbind(rep2TH0hpitotal_exons,rep2TH4hpitotal_exons,rep2TH8hpitotal_exons)
  colnames(rep2total_exons1) <-c("rep2total_exonsOnehour","rep2total_exonsFourhour","rep2total_exonsEighthour")
  
  rep2TH0hpitotal_introns=rep2H0hpiRPKMsOut$rpkms$total_introns
  rep2TH4hpitotal_introns=rep2H4hpiRPKMsOut$rpkms$total_introns
  rep2TH8hpitotal_introns=rep2H8hpiRPKMsOut$rpkms$total_introns
  rep2total_introns1=cbind(rep2TH0hpitotal_introns,rep2TH4hpitotal_introns,rep2TH8hpitotal_introns)
  colnames(rep2total_introns1) <-c("rep2total_intronsOnehour","rep2total_intronsFourhour","rep2total_intronsEighthour")
  
  tpts <- c(0,4,8,0,4,8)
  tL <- 4
  foursu_exonsall=cbind(rep1foursu_exons,rep2foursu_exons1)
  colnames(foursu_exonsall) <-c("foursu_exonsREP1Onehour","foursu_exonsREP1Fourhour","foursu_exonsREP1Eighthour","foursu_exonsREP2Onehour","foursu_exonsREP2Fourhour","foursu_exonsREP2Eighthour")
  
  total_exonsall=cbind(rep1total_exons,rep2total_exons1)
  colnames(total_exonsall) <-c("total_exonsREP1Onehour","total_exonsREP1Fourhour","total_exonsREP1Eighthour","total_exonsREP2Onehour","total_exonsREP2Fourhour","total_exonsREP2Eighthour")
  
  foursu_intronsall=cbind(rep1foursu_introns,rep2foursu_introns1)
  colnames(foursu_intronsall) <-c("foursu_intronsREP1Onehour","foursu_intronsREP1Fourhour","foursu_intronsREP1Eighthour","foursu_intronsREP2Onehour","foursu_intronsREP2Fourhour","foursu_intronsREP2Eighthour")
  
  total_intronsall=cbind(rep1total_introns,rep2total_introns1)
  colnames(total_intronsall) <-c("total_intronsREP1Onehour","total_intronsREP1Fourhour","total_intronsREP1Eighthour","total_intronsREP2Onehour","total_intronsREP2Fourhour","total_intronsREP2Eighthour")
  
  # GeneLabAREP12<-newINSPEcT(tpts,tL,foursu_exonsall, total_exons,foursu_intronsall,total_intronsall,BPPARAM=SerialParam())

  result=newINSPEcT(tpts, tL, foursu_exonsall, total_exonsall,
                                              foursu_intronsall, total_intronsall, BPPARAM=SerialParam())
  Express=result@ratesFirstGuess@assayData$exprs
  return(Express[,7:length(Express[1,])])
}


# get the count result
GetCountNumberResult<-function()
{
library(GenomicFeatures)
library(INSPEcT)
library(Rsamtools)
library(deSolve)
library(BiocParallel)
library(compiler)
library(sqldf)
txdb3 <- makeTxDbFromGFF(file="D:/ivan/Harm/NewTestOne/NewTest/Newgencode.v23.PR8.gtf",format="gtf",taxonomyId=9606)
# load("D:/ivan/Harm/NewTestFive/tran/InputData/count.RData")
# delvari("^C")
CountLabArep1=GetcountResultInspect("LabA0hpi.bam",'UnlA0hpi.bam','LabA4hpi.bam','UnlA4hpi.bam','LabA8hpi.bam','UnlA8hpi.bam',"Rep1","gene")
colnames(CountLabArep1)=c("id","LabA0","LabA4","LabA8","LabA10","LabA11","LabA12","UnlA0","UnlA4","UnlA8","UnlA10","UnlA11","UnlA12")
CountLabArep2=GetcountResultInspect("LabA0.bam",'ULabA0.bam','LabA4.bam','ULabA4.bam','LabA8.bam','ULabA8.bam',"Rep2","gene")
colnames(CountLabArep2)=c("id","LabA0_2","LabA4_2","LabA8_2","LabA10_2","LabA11_2","LabA12_2","UnlA0_2","UnlA4_2","UnlA8_2","UnlA10_2","UnlA11_2","UnlA12_2")

CountLabCrep1=GetcountResultInspect("LabC0hpi.bam",'UnlC0hpi.bam','LabC4hpi.bam','UnlC4hpi.bam','LabC8hpi.bam','UnlC8hpi.bam',"Rep1","gene")
colnames(CountLabCrep1)=c("id","LabC0","LabC4","LabC8","LabC10","LabC11","LabC12","UnlC0","UnlC4","UnlC8","UnlC10","UnlC11","UnlC12")
CountLabCrep2=GetcountResultInspect("LabC0.bam",'ULabC0.bam','LabC4.bam','ULabC4.bam','LabC8.bam','ULabC8.bam',"Rep2","gene")
colnames(CountLabCrep2)=c("id","LabC0_2","LabC4_2","LabC8_2","LabC10_2","LabC11_2","LabC12_2","UnlC0_2","UnlC4_2","UnlC8_2","UnlC10_2","UnlC11_2","UnlC12_2")

restult=merge(CountLabArep1, CountLabArep2, by ="id",all = TRUE,incomparables = NULL)
restult=merge(restult, CountLabCrep1, by ="id",all = TRUE,incomparables = NULL) 
restult=merge(restult, CountLabCrep2, by ="id",all = TRUE,incomparables = NULL) 

return(restult)
}

