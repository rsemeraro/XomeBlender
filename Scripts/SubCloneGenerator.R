vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,","))

##  Setting input paths and parameters ###
FileVCFIn <- as.character(split.vars[1])
Label <- as.character(split.vars[2])
Model <- as.character(split.vars[3])
TotalMutation <- as.numeric(split.vars[4])
NClone <- as.numeric(split.vars[5])
Path2VCF <- as.character(split.vars[6])

TableVCF <- read.table(FileVCFIn, header = F, sep ="\t",quote="",fill=TRUE)
ChrVCF<-as.character(TableVCF[,1])
PosVCF<-as.numeric(TableVCF[,2])
RefVCF<-as.character(TableVCF[,4])
AltVCF<-as.character(TableVCF[,5])

StringFormat<-'##fileformat=VCFv4.0'
StringDate<-paste("##fileDate=",format(Sys.time(), "%Y%d%m"),sep="")
Header<-c(StringFormat,StringDate)
HeaderMat<-c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT')


######################
#### Linear Model ####
######################
if (Model=="Linear")
{
  indAll<-c(1:nrow(TableVCF))
  indControl<-sort(sample(indAll,TotalMutation))
  VCF2Remove<-TableVCF[indControl,]
  VCF2VariantClone<-TableVCF[indControl,]
  FileVCF2Remove<-file.path(Path2VCF,paste(Label,"_Control.remove",sep=""))
  FileVCFVariants<-file.path(Path2VCF,paste(Label,"_Subclone_1.variants",sep=""))
  write.table(VCF2Remove,FileVCF2Remove , sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
  write.table(VCF2VariantClone, FileVCFVariants, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
  
  #### Saving Hetero Variants VCF for CNV creation ####
  ExpLabelOut<-paste(Label,"_Subclone_1",sep="")
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
      
      FileVCF2Remove<-file.path(Path2VCF,paste(Label,"_Subclone_",i,".remove",sep=""))
      FileVCFVariants<-file.path(Path2VCF,paste(Label,"_Subclone_",i,".variants",sep=""))
      
      write.table(VCF2Remove,FileVCF2Remove , sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
      write.table(VCF2VariantClone, FileVCFVariants, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
      
      #### Saving Hetero Variants VCF for CNV creation ####
      ExpLabelOut<-paste(Label,"_Subclone_",i,sep="")
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
    }
  }
}



######################
#### Branched Model ####
######################
if (Model=="Branched")
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
  
  
  FileVCF2Remove<-file.path(Path2VCF,paste(Label,"_Control.remove",sep=""))
  FileVCF2Remove1<-file.path(Path2VCF,paste(Label,"_Subclone_1.remove",sep=""))
  FileVCF2Remove2<-file.path(Path2VCF,paste(Label,"_Subclone_2.remove",sep=""))
  FileVCFVariants1<-file.path(Path2VCF,paste(Label,"_Subclone_1.variants",sep=""))
  FileVCFVariants2<-file.path(Path2VCF,paste(Label,"_Subclone_2.variants",sep=""))
  
  write.table(VCF2Remove,FileVCF2Remove , sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
  write.table(VCF2Remove1,FileVCF2Remove1 , sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
  write.table(VCF2Remove2,FileVCF2Remove2 , sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
  
  write.table(VCF2VariantClone1, FileVCFVariants1, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
  write.table(VCF2VariantClone2, FileVCFVariants2, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
  
  #### Saving Hetero Variants VCF for CNV creation from SublClone1 ####
  TableVCF1<-TableVCF[-indSubClone1,]
  ExpLabelOut1<-paste(Label,"_Subclone_1",sep="")
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
  
  #### Saving Hetero Variants VCF for CNV creation from SublClone2 ####
  TableVCF2<-TableVCF[-indSubClone2,]
  ExpLabelOut2<-paste(Label,"_Subclone_2",sep="")
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
  
  
  if (NClone==3)
  {
    
    
    NSplit<-round((length(indControl)/4))
    
    
    indSubClone3<-indSubClone1[1:NSplit]
    indSubClone4<-indSubClone1[(NSplit+1):length(indSubClone1)]
    
    VCF2Remove3<-TableVCF[indSubClone4,]
    
    
    VCF2VariantClone3<-TableVCF[indSubClone3,]
    
    
    
    
    FileVCF2Remove3<-file.path(Path2VCF,paste(Label,"_Subclone_3",".remove",sep=""))
    FileVCFVariants3<-file.path(Path2VCF,paste(Label,"_Subclone_3",".variants",sep=""))
    
    write.table(VCF2Remove3,FileVCF2Remove3 , sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
    write.table(VCF2VariantClone3, FileVCFVariants3, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
    
    #### Saving Hetero Variants VCF for CNV creation ####
    TableVCF3<-TableVCF1[-indSubClone3,]
    ExpLabelOut1<-paste(Label,"_Subclone_3",sep="")
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
    
    
    
    FileVCF2Remove3<-file.path(Path2VCF,paste(Label,"_Subclone_3",".remove",sep=""))
    FileVCFVariants3<-file.path(Path2VCF,paste(Label,"_Subclone_3",".variants",sep=""))
    FileVCF2Remove4<-file.path(Path2VCF,paste(Label,"_Subclone_4",".remove",sep=""))
    FileVCFVariants4<-file.path(Path2VCF,paste(Label,"_Subclone_4",".variants",sep=""))
    
    
    write.table(VCF2Remove3,FileVCF2Remove3 , sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
    write.table(VCF2VariantClone3, FileVCFVariants3, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
    write.table(VCF2Remove4,FileVCF2Remove4 , sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
    write.table(VCF2VariantClone4, FileVCFVariants4, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
    
    #### Saving Hetero Variants VCF for CNV creation ####
    TableVCF3<-TableVCF[-indSubClone3,]
    ExpLabelOut1<-paste(Label,"_Subclone_3",sep="")
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
    
    
    #### Saving Hetero Variants VCF for CNV creation ####
    TableVCF4<-TableVCF[-indSubClone5,]
    ExpLabelOut1<-paste(Label,"_Subclone_4",sep="")
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
    
  }
}
