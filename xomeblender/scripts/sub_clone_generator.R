vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,","))
FinalRData <- split.vars[1]

load(FinalRData)
TableVCF <- read.table(FileVCFIn, header = F, sep ="\t",quote="",fill=TRUE)
ChrVCF<-as.character(TableVCF[,1])
PosVCF<-as.numeric(TableVCF[,2])
RefVCF<-as.character(TableVCF[,4])
AltVCF<-as.character(TableVCF[,5])


StringHeader<-paste("cat ",FileVCFIn," | head -200 | grep '#'",sep="")
TotalHeader<-system(StringHeader,intern = TRUE)
NLine<-length(TotalHeader)
Header<-TotalHeader[-NLine]
HeaderMat<-c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT')


######################
#### Linear Model ####
######################

if (NClone < 2)
{
  Model <- "linear"
}

if (Model=="linear")
{
  indAll<-c(1:nrow(TableVCF))
  indControl<-sort(sample(indAll,TotalMutation))
  VCF2Remove<-TableVCF[indControl,]
  VCF2VariantClone<-TableVCF[indControl,]
  FileVCF2Remove<-file.path(Path2VCF,paste(Label,"_C.remove",sep=""))
  FileVCFVariants<-file.path(Path2VCF,paste(Label,"_S1.variants",sep=""))
  write.table(VCF2Remove,FileVCF2Remove , sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
  write(paste(Label,"_C.remove",sep=""), stdout())
  write.table(VCF2VariantClone, FileVCFVariants, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
  
  #### Saving Hetero Variants VCF for CNV creation ####
  ExpLabelOut<-paste(Label,"_S1",sep="")
  HeadMatSample<-c(HeaderMat,ExpLabelOut)
  TableVCF2Save<-c()
  for (zz in 1:10)
  {
    TableVCF2Save<-cbind(TableVCF2Save,as.character(TableVCF[,zz]))
    
  }
  MatOut<-rbind(HeadMatSample,TableVCF2Save)
  
  FileOut<-file.path(Path2VCF,paste(ExpLabelOut,".vcf",sep=""))
  
  zz <- file(FileOut, "w")
  cat(Header,file=zz,sep = "\n")
  write(t(MatOut),file=zz,sep = "\t",ncolumns=ncol(MatOut),append=TRUE)
  close(zz)
  compressedout <- paste(ExpLabelOut,".vcf.gz",sep="")
  BgzipString <- paste("bgzip -c", FileOut, ">", compressedout, sep=" ")
  system(BgzipString, wait = TRUE)
  IndexingString <- paste("tabix -p vcf", compressedout, sep=" ")
  system(IndexingString, wait = TRUE)
  
  if (NClone>1)
  {
    Mutation2Sub<-round(TotalMutation/NClone)
    Mutation2Keep<-TotalMutation-Mutation2Sub
    Spick<-1
    Epick<-Mutation2Sub
    for (i in 2:NClone)
    {
      indPickSub<-c(Spick:(Epick*(i-1)))
      indPickVar<-c((Epick*(i-1)+1):TotalMutation)
      indSub<-indControl[indPickSub]
      indVar<-indControl[indPickVar]
      VCF2Remove<-TableVCF[indSub,]
      VCF2VariantClone<-TableVCF[indVar,]
      
      FileVCF2Remove<-file.path(Path2VCF,paste(Label,"_S",i,".remove",sep=""))
      FileVCFVariants<-file.path(Path2VCF,paste(Label,"_S",i,".variants",sep=""))
      
      write.table(VCF2Remove,FileVCF2Remove , sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
      write(paste(Label,"_S",i,".remove",sep=""), stdout())
      write.table(VCF2VariantClone, FileVCFVariants, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
      
      #### Saving Hetero Variants VCF for CNV creation ####
      ExpLabelOut<-paste(Label,"_S",i,sep="")
      HeadMatSample<-c(HeaderMat,ExpLabelOut)
      TableVCF2Save<-c()
      for (zz in 1:10)
      {
        TableVCF2Save<-cbind(TableVCF2Save,as.character(TableVCF[-sort(indSub),zz]))
        
      }
      MatOut<-rbind(HeadMatSample,TableVCF2Save)
      
      FileOut<-file.path(Path2VCF,paste(ExpLabelOut,".vcf",sep=""))
      
      zz <- file(FileOut, "w")
      cat(Header,file=zz,sep = "\n")
      write(t(MatOut),file=zz,sep = "\t",ncolumns=ncol(MatOut),append=TRUE)
      close(zz)
      compressedout <- paste(ExpLabelOut,".vcf.gz",sep="")
      BgzipString <- paste("bgzip -c", FileOut, ">", compressedout, sep=" ")
      system(BgzipString, wait = TRUE)
      IndexingString <- paste("tabix -p vcf", compressedout, sep=" ")
      system(IndexingString, wait = TRUE)
    }
  }
}



######################
#### Branched Model ####
######################
if (Model=="branched")
{  
  indAll<-c(1:nrow(TableVCF))
  indControl<-sort(sample(indAll,TotalMutation))
  NSplit<-round((length(indControl)/2))
  VCF2Remove<-TableVCF[indControl,]
  
  
  indSubClone1<-indControl[1:NSplit]
  indSubClone2<-indControl[(NSplit+1):length(indControl)]
  
  
  VCF2Remove1<-TableVCF[indSubClone2,]
  VCF2Remove2<-TableVCF[indSubClone1,]
  
  VCF2VariantClone1<-TableVCF[indSubClone1,]
  VCF2VariantClone2<-TableVCF[indSubClone2,]
  
  
  FileVCF2Remove<-file.path(Path2VCF,paste(Label,"_C.remove",sep=""))
  write(paste(Label,"_C.remove",sep=""), stdout())
  FileVCF2Remove1<-file.path(Path2VCF,paste(Label,"_S1.remove",sep=""))
  write(paste(Label,"_S1.remove",sep=""), stdout())
  FileVCF2Remove2<-file.path(Path2VCF,paste(Label,"_S2.remove",sep=""))
  write(paste(Label,"_S2.remove",sep=""), stdout())
  FileVCFVariants1<-file.path(Path2VCF,paste(Label,"_S1.variants",sep=""))
  FileVCFVariants2<-file.path(Path2VCF,paste(Label,"_S2.variants",sep=""))
  
  write.table(VCF2Remove,FileVCF2Remove , sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
  write.table(VCF2Remove1,FileVCF2Remove1 , sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
  write.table(VCF2Remove2,FileVCF2Remove2 , sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
  
  write.table(VCF2VariantClone1, FileVCFVariants1, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
  write.table(VCF2VariantClone2, FileVCFVariants2, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
  
  #### Saving Hetero Variants VCF for CNV creation from SublClone1 ####
  TableVCF1<-TableVCF[-indSubClone1,]
  ExpLabelOut1<-paste(Label,"_S1",sep="")
  HeadMatSample<-c(HeaderMat,ExpLabelOut1)
  TableVCF2Save<-c()
  for (zz in 1:10)
  {
    TableVCF2Save<-cbind(TableVCF2Save,as.character(TableVCF1[,zz]))
    
  }
  MatOut<-rbind(HeadMatSample,TableVCF2Save)
  
  FileOut1<-file.path(Path2VCF,paste(ExpLabelOut1,".vcf",sep=""))
  
  zz <- file(FileOut1, "w")
  cat(Header,file=zz,sep = "\n")
  write(t(MatOut),file=zz,sep = "\t",ncolumns=ncol(MatOut),append=TRUE)
  close(zz)
  compressedout1 <- paste(ExpLabelOut1,".vcf.gz",sep="")
  BgzipString <- paste("bgzip -c", FileOut1, ">", compressedout1, sep=" ")
  system(BgzipString, wait = TRUE)
  IndexingString <- paste("tabix -p vcf", compressedout1, sep=" ")
  system(IndexingString, wait = TRUE)
  #### Saving Hetero Variants VCF for CNV creation from SublClone2 ####
  TableVCF2<-TableVCF[-indSubClone2,]
  ExpLabelOut2<-paste(Label,"_S2",sep="")
  HeadMatSample<-c(HeaderMat,ExpLabelOut2)
  TableVCF2Save<-c()
  for (zz in 1:10)
  {
    TableVCF2Save<-cbind(TableVCF2Save,as.character(TableVCF2[,zz]))
    
  }
  MatOut<-rbind(HeadMatSample,TableVCF2Save)
  
  FileOut2<-file.path(Path2VCF,paste(ExpLabelOut2,".vcf",sep=""))
  
  zz <- file(FileOut2, "w")
  cat(Header,file=zz,sep = "\n")
  write(t(MatOut),file=zz,sep = "\t",ncolumns=ncol(MatOut),append=TRUE)
  close(zz)
  compressedout2 <- paste(ExpLabelOut2,".vcf.gz",sep="")
  BgzipString <- paste("bgzip -c", FileOut2, ">", compressedout2, sep=" ")
  system(BgzipString, wait = TRUE)
  IndexingString <- paste("tabix -p vcf", compressedout2, sep=" ")
  system(IndexingString, wait = TRUE)
  
  if (NClone==3)
  {
    
    
    NSplit<-round((length(indControl)/4))
    
    
    indSubClone3<-indSubClone1[1:NSplit]
    indSubClone4<-indSubClone1[(NSplit+1):length(indSubClone1)]
    
    VCF2Remove3<-TableVCF[indSubClone4,]
    
    
    VCF2VariantClone3<-TableVCF[indSubClone3,]
    
    
    
    
    FileVCF2Remove3<-file.path(Path2VCF,paste(Label,"_S3",".remove",sep=""))
    FileVCFVariants3<-file.path(Path2VCF,paste(Label,"_S3",".variants",sep=""))
    write(paste(Label,"_S3.remove",sep=""), stdout())
    write.table(VCF2Remove3,FileVCF2Remove3 , sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
    write.table(VCF2VariantClone3, FileVCFVariants3, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
    
    #### Saving Hetero Variants VCF for CNV creation ####
    TableVCF3<-TableVCF1[-indSubClone3,]
    ExpLabelOut1<-paste(Label,"_S3",sep="")
    HeadMatSample<-c(HeaderMat,ExpLabelOut1)
    TableVCF2Save<-c()
    for (zz in 1:10)
    {
      TableVCF2Save<-cbind(TableVCF2Save,as.character(TableVCF3[,zz]))
      
    }
    MatOut<-rbind(HeadMatSample,TableVCF2Save)
    
    FileOut1<-file.path(Path2VCF,paste(ExpLabelOut1,".vcf",sep=""))
    
    zz <- file(FileOut1, "w")
    cat(Header,file=zz,sep = "\n")
    write(t(MatOut),file=zz,sep = "\t",ncolumns=ncol(MatOut),append=TRUE)
    close(zz)
    compressedout1 <- paste(ExpLabelOut1,".vcf.gz",sep="")
    BgzipString <- paste("bgzip -c", FileOut1, ">", compressedout1, sep=" ")
    system(BgzipString, wait = TRUE)
    IndexingString <- paste("tabix -p vcf", compressedout1, sep=" ")
    system(IndexingString, wait = TRUE)
  }
  
  
  if (NClone==4)
  {
    
    
    NSplit<-round((length(indControl)/4))
    indSubClone5<-indSubClone2[1:NSplit]
    indSubClone6<-indSubClone2[(NSplit+1):length(indSubClone2)]
    
    indSubClone3<-indSubClone1[1:NSplit]
    indSubClone4<-indSubClone1[(NSplit+1):length(indSubClone1)]
    
    VCF2Remove3<-TableVCF[indSubClone4,]
    VCF2Remove4<-TableVCF[indSubClone6,]
    
    VCF2VariantClone3<-TableVCF[indSubClone3,]
    VCF2VariantClone4<-TableVCF[indSubClone5,]
    
    
    
    FileVCF2Remove3<-file.path(Path2VCF,paste(Label,"_S3",".remove",sep=""))
    write(paste(Label,"_S3.remove",sep=""), stdout())
    FileVCFVariants3<-file.path(Path2VCF,paste(Label,"_S3",".variants",sep=""))
    FileVCF2Remove4<-file.path(Path2VCF,paste(Label,"_S4",".remove",sep=""))
    write(paste(Label,"_S.remove",sep=""), stdout())
    FileVCFVariants4<-file.path(Path2VCF,paste(Label,"_S4",".variants",sep=""))
    
    
    write.table(VCF2Remove3,FileVCF2Remove3 , sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
    write.table(VCF2VariantClone3, FileVCFVariants3, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
    write.table(VCF2Remove4,FileVCF2Remove4 , sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
    write.table(VCF2VariantClone4, FileVCFVariants4, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
    
    #### Saving Hetero Variants VCF for CNV creation ####
    TableVCF3<-TableVCF[-indSubClone3,]
    ExpLabelOut1<-paste(Label,"_S3",sep="")
    HeadMatSample<-c(HeaderMat,ExpLabelOut1)
    TableVCF2Save<-c()
    for (zz in 1:10)
    {
      TableVCF2Save<-cbind(TableVCF2Save,as.character(TableVCF3[,zz]))
      
    }
    MatOut<-rbind(HeadMatSample,TableVCF2Save)
    
    FileOut1<-file.path(Path2VCF,paste(ExpLabelOut1,".vcf",sep=""))
    
    zz <- file(FileOut1, "w")
    cat(Header,file=zz,sep = "\n")
    write(t(MatOut),file=zz,sep = "\t",ncolumns=ncol(MatOut),append=TRUE)
    close(zz)
    compressedout1 <- paste(ExpLabelOut1,".vcf.gz",sep="")
    BgzipString <- paste("bgzip -c", FileOut1, ">", compressedout1, sep=" ")
    system(BgzipString, wait = TRUE)
    IndexingString <- paste("tabix -p vcf", compressedout1, sep=" ")
    system(IndexingString, wait = TRUE)
    
    #### Saving Hetero Variants VCF for CNV creation ####
    TableVCF4<-TableVCF[-indSubClone5,]
    ExpLabelOut1<-paste(Label,"_S4",sep="")
    HeadMatSample<-c(HeaderMat,ExpLabelOut1)
    TableVCF2Save<-c()
    for (zz in 1:10)
    {
      TableVCF2Save<-cbind(TableVCF2Save,as.character(TableVCF4[,zz]))
      
    }
    MatOut<-rbind(HeadMatSample,TableVCF2Save)
    
    FileOut1<-file.path(Path2VCF,paste(ExpLabelOut1,".vcf",sep=""))
    
    zz <- file(FileOut1, "w")
    cat(Header,file=zz,sep = "\n")
    write(t(MatOut),file=zz,sep = "\t",ncolumns=ncol(MatOut),append=TRUE)
    close(zz)
    compressedout1 <- paste(ExpLabelOut1,".vcf.gz",sep="")
    BgzipString <- paste("bgzip -c", FileOut1, ">", compressedout1, sep=" ")
    system(BgzipString, wait = TRUE)
    IndexingString <- paste("tabix -p vcf", compressedout1, sep=" ")
    system(IndexingString, wait = TRUE)
  }
}
