# inputfile="D:/guojun/1_R1_001.fastq"
# seqone="CTGGAATTCGCGGTTAAA"
# # seqtwo="AGTAGAAACAAGGTCGTTTTT"
# seqtwo="CTGGAATTCGCGGTTAAA"
# outfile="D:/guojun/result.fastq"
Esubset<-function(inputfile,outfile,seqone,seqtwo)
{
  text <- readLines(inputfile, encoding = "UTF-8")
  firstrecord=as.data.frame(text[seq(1,length(text),4)])
  secondrecord=as.data.frame(text[seq(2,length(text),4)])
  keptrecord=cbind(firstrecord,secondrecord)
  index1=as.data.frame(grep(seqone, as.character(keptrecore[,2])))
  index2=as.data.frame(grep(seqtwo, as.character(keptrecore[,2])))
  index=as.numeric(unlist(intersect(index1,index2))) 
  SearchRecord=keptrecore[index,]
  write.table(SearchRecord, file =outfile, append = FALSE, quote = F, sep = "\t", row.names =F,col.names = F)

  
}

Esubset("D:/guojun/1_R1_001.fastq","D:/guojun/result.fastq","CTGGAATTCGCGGTTAAA","CTGGAATTCGCGGTTAAA")
