### Open input ###

vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars, ","))
PathToInput <- split.vars[1]
OriginalVcf <- read.table(PathToInput,sep="\t",header=F,quote = "\"",fill=T)
CovThr <- 10
HomoRefThr <- 0.01
EtheroAltThr <- 0.1

### Header Parsing ###

StringHeader<-paste(" head -200 ",PathToInput," | grep '#'",sep="")
TotalHeader<-system(StringHeader,intern = TRUE)
NLine<-length(TotalHeader)
HeaderNames<-TotalHeader[NLine]
Header<-TotalHeader[-NLine]
HeaderSplit<-unlist(strsplit(HeaderNames,"\t"))
ColumnVec <- c(10:length(HeaderSplit))

### Alternative aplotypes removing ###

AltAllele<-as.character(OriginalVcf[,5])
indMulti<-grep(",",AltAllele)
if (length(indMulti)!=0)
{
  NoAltAploVcf <- OriginalVcf[-indMulti,]
}

if (length(indMulti)==0)
{
  NoAltAploVcf <- OriginalVcf
}

### Filtering ###

NoAltAploVcfOut <- c()

for (c in 1:length(ColumnVec))
{
### Coverage Filtering ###

  NoAltAploVcfNew<-as.character(NoAltAploVcf[,(ColumnVec[c])])
  Fixer <- which(substr(NoAltAploVcfNew,1 ,3) == "./.")
  NoAltAploVcfNew[Fixer] <- "./."
  GoodPos <- which(substr(NoAltAploVcfNew,1 ,3) != "./.")
  AllGenoVec <- NoAltAploVcfNew[GoodPos]
  SplitVec <- unlist(strsplit(AllGenoVec, ":"))
  NN <- length(SplitVec)
  indAD <- seq(2, NN, by = 5)
  AllelicDepths <- as.character(SplitVec[indAD])
  SplitVec2 <- unlist(strsplit(AllelicDepths, ","))
  indRef <- seq(1, length(SplitVec2), by=2)
  indAlt <- seq(2, length(SplitVec2), by=2)
  Ref <- as.numeric(SplitVec2[indRef])
  Alt <- as.numeric(SplitVec2[indAlt])
  Cov <- Ref+Alt
  CovFilter <- which(Cov < CovThr)
  NoAltAploVcfNew[GoodPos[CovFilter]]<-"./."

  ### Bad 0/0 removing ###

  Bait <- which(substr(NoAltAploVcfNew,1 ,3) == "0/0")
  ZeroZeroVec <- as.character(NoAltAploVcfNew[Bait])
  SplitVec <- unlist(strsplit(ZeroZeroVec, ":"))
  NN <- length(SplitVec)
  indAD <- seq(2, NN, by = 5)
  AllelicDepths <- as.character(SplitVec[indAD])
  SplitVec2 <- unlist(strsplit(AllelicDepths, ","))
  indRef <- seq(1, length(SplitVec2), by=2)
  indAlt <- seq(2, length(SplitVec2), by=2)
  Ref <- as.numeric(SplitVec2[indRef])
  Alt <- as.numeric(SplitVec2[indAlt])
  BAF <- (Alt/(Alt+Ref))
  RatioEvaluation <- which(as.numeric(BAF) >= HomoRefThr)
  NoAltAploVcfNew[Bait[RatioEvaluation]]<-"./."

  ### Bad 0/1 and 1/1 removing ###

  Bait <- which(substr(NoAltAploVcfNew,1 ,3) == "0/1" | substr(NoAltAploVcfNew,1 ,3) == "1/1")
  EtheroAltVec <- as.character(NoAltAploVcfNew[Bait])
  SplitVec <- unlist(strsplit(EtheroAltVec, ":"))
  NN <- length(SplitVec)
  indAD <- seq(2, NN, by = 5)
  AllelicDepths <- as.character(SplitVec[indAD])
  SplitVec2 <- unlist(strsplit(AllelicDepths, ","))
  indRef <- seq(1, length(SplitVec2), by=2)
  indAlt <- seq(2, length(SplitVec2), by=2)
  Ref <- as.numeric(SplitVec2[indRef])
  Alt <- as.numeric(SplitVec2[indAlt])
  BAF <- (Alt/(Alt+Ref))
  RatioEvaluation <- which(as.numeric(BAF) < EtheroAltThr)
  NoAltAploVcfNew[Bait[RatioEvaluation]]<-"./."

  NoAltAploVcfOut<-cbind(NoAltAploVcfOut,NoAltAploVcfNew)
}

NoAltAploVcfSave<-cbind(NoAltAploVcf[,c(1:9)],NoAltAploVcfOut)

### Saving File ###

cat(TotalHeader,file=PathToInput,sep = "\n")
write(t(NoAltAploVcfSave),file=PathToInput,sep = "\t",ncolumns=ncol(NoAltAploVcfSave),append=TRUE)
